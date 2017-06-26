#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>
#include "boost/multi_array.hpp"
#include "saver.hpp"
#include <vexcl/vexcl.hpp>
#include "pr_util.h"
#include "omp.h"
#include <cassert>

#ifdef nDEBUG
#  undef nDEBUG
#endif

/*AMGCL include section*/
#include <amgcl/amg.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/backend/vexcl_static_matrix.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/io/mm.hpp>
/*end of section*/


namespace amgcl
{
    profiler<> prof("project suffosion");
}
using amgcl::prof;

int main(int argc, char *argv[])
{
    typedef amgcl::static_matrix<double, 3, 3> mat_type;
    typedef amgcl::static_matrix<double, 3, 1> vec_type;

    // Initialize VexCL context.
    vex::Context ctx( vex::Filter::Env && vex::Filter::DoublePrecision );

    if (!ctx)
    {
        std::cerr << "no GPUs" << std::endl;
        return 1;
    }

    std::cout << ctx << std::endl;

    vex::scoped_program_header header(ctx,
            amgcl::backend::vexcl_static_matrix_declaration<double,3>());

    //Vectors of variables
    std::vector<double> pressure(n, p_top);
    std::vector<cell> flow(n); //velocities, x,y,z components
    std::vector<cell> c_vol(n); //volume concentration
    std::vector<double> source(n, 0.);//source of solid phase
    std::vector<double> concentration(n, c_0); 
    std::vector<double> tmp_conc(n, c_0); 
    std::vector<double> permeability(n, k_0);
    std::vector<double> matrix(n, 0.); //matrix of permeability
    std::vector<double> u_x(n, 0.), u_y(n, 0.), u_z(n, 0.); //displacements
    std::vector<double> dilatation(n, 0.); 
    std::vector<double> dil_dt(n, 0.); //dilatation time derivative
    std::vector<double> porosity(n, fi_0);
    std::vector<double> por_dt(n, 0.); //porosity time derivative
    std::vector<int> wells; //coordinates of wells

    //Vectors for calculations - pressure
    std::vector<double> val_pr; //values of nonzero entries
    std::vector<int>    col_pr; //column numbers of nonzero entries
    std::vector<int>    ptr_pr; //points to the start of each row in col_pr
    std::vector<double> rhs_pr(n, 0.); //right-hand side 
    vex::vector<double> rhs_dev_pr(ctx.queue(), rhs_pr); //rhs on device
    vex::vector<double> x_pr(ctx.queue(), pressure); //unknown vector on device

    //Vectors for calculations - displacements
    std::vector<double> val_u; //values of nonzero entries
    std::vector<int>    col_u; //column numbers of nonzero entries
    std::vector<int>    ptr_u; //points to the start of each row in col_pr 
    std::vector<double> rhs_u(3 * n, 0.); // right-hand side 
    std::vector<double> disp(3 * n, 0.); 
    vex::vector<vec_type> rhs_dev_u(ctx.queue(), n); // rhs on device
    vex::vector<vec_type> x_u(ctx.queue(), n);
    amgcl::backend::clear(x_u);
    std::string filename;
    
    std::vector<double> temp_val(n, 0.); 

    std::map<std::string, double*> save_data = 
    {
        {"pressure", pressure.data()},
        {"dilatation", dilatation.data()},
        {"porosity", porosity.data()},
        {"temp", temp_val.data()},
        {"permeability", permeability.data()},
        {"dil_dt", dil_dt.data()},
        {"concentration", concentration.data()},
        {"source", source.data()},
        {"u_x", u_x.data()},
        {"u_y", u_y.data()},
        {"u_z", u_z.data()}
    };
    saver data_saver("suffosion", n_x, n_y, n_z, h);

    //other variables
    const auto duration = ( (argc > 1) ? std::stoul( argv[1] ) : 1 ); 
    int writer_step = 99;

    add_well(n_x/2, n_y/2, wells);

    std::sort(wells.begin(), wells.end());

    // Define the AMG type:
    typedef amgcl::backend::vexcl<double>   SBackend;
    typedef amgcl::backend::vexcl<mat_type> BBackend;

    typedef amgcl::make_solver<
        amgcl::amg<
        BBackend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::damped_jacobi
            >,
        amgcl::solver::bicgstab< BBackend >
            > USolver;

    // Define the AMG type:
    typedef amgcl::make_solver<
        amgcl::amg<
        SBackend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
            >,
        amgcl::solver::gmres< SBackend >
            > PSolver;

    col_pr.reserve(7 * n);
    val_pr.reserve(7 * n);
    ptr_pr.reserve(n);

    int    iters;
    double error;

    prof.tic("build press matrix");
    build_press_mat(col_pr, val_pr, ptr_pr, rhs_pr, permeability, wells,
            source, dil_dt);

    PSolver solve_press( boost::tie(n, ptr_pr, col_pr, val_pr) );
    prof.toc("build press matrix");

    prof.tic("solve press matrix");
    vex::copy(rhs_pr, rhs_dev_pr);
    boost::tie(iters, error) = solve_press(rhs_dev_pr, x_pr);
    vex::copy(x_pr, pressure);   
    prof.toc("solve press matrix");

    std::cout << "iters P = " << iters << std::endl;
    std::cout << "error P = " << error << std::endl;

    int n_u = 3 * n;
    col_u.reserve(5 * n_u);
    val_u.reserve(5 * n_u);
    ptr_u.reserve(n_u);

    prof.tic("build disp matrix");
    build_disp_mat(col_u, val_u, ptr_u, wells);
    fill_disp_rhs(pressure, rhs_u, wells);
    USolver solve_disp( amgcl::adapter::block_matrix<3, mat_type>(
                boost::tie(n_u, ptr_u, col_u, val_u) ) );
    prof.toc("build disp matrix");

    prof.tic("solve disp");

    vec_type const * fptr = reinterpret_cast<vec_type const *>(&rhs_u[0]);
    vec_type       * xptr = reinterpret_cast<vec_type       *>(&disp[0]);

    vex::copy(fptr, fptr + n, rhs_dev_u.begin());
    boost::tie(iters, error) = solve_disp(rhs_dev_u, x_u);
    vex::copy(x_u.begin(), x_u.end(), xptr);
    prof.toc("solve disp");

    std::cout << "iters U = " << iters << std::endl;
    std::cout << "error U = " << error << std::endl;

    dil_calc(disp, dilatation, dil_dt, wells);
    std::fill(dil_dt.begin(), dil_dt.end(), 0.);

    prof.tic("time cycle");
    for(auto t = 0; t < 0; ++t)
    {
        //pressure
        build_press_mat(col_pr, val_pr, ptr_pr, rhs_pr, permeability, wells,
                source, dil_dt);
        vex::copy(rhs_pr, rhs_dev_pr);
        PSolver solve_press(boost::tie(n, ptr_pr, col_pr, val_pr));
        solve_press(rhs_dev_pr, x_pr);
        vex::copy(x_pr, pressure);

        //displacements
        fill_disp_rhs(pressure, rhs_u, wells);
        vex::copy(fptr, fptr + n, rhs_dev_u.begin());
        solve_disp(rhs_dev_u, x_u);
        vex::copy(x_u.begin(), x_u.end(), xptr);
       
        for(auto index = 0; index < n; ++index)
        {
            u_x[index] = disp[3*index + 0];
            u_y[index] = disp[3*index + 1];
            u_z[index] = disp[3*index + 2];
        }

        //dilatation
        dil_calc(disp, dilatation, dil_dt, wells);

        //velocity

        //concentration

        //porosity&source

        //permeability

        //write results
        if( (writer_step++) == 100 || t == duration - 1 )
        {
            data_saver.add_step(t*h_t, save_data);
            
            writer_step = 0;
        }
    }

    for(auto index = 0; index < n; ++index)
        if(is_well(index, wells))
            concentration[index] = 0.1;

    get_flow(flow, pressure, permeability);

    for(auto t = 0; t < 100; ++t)
        lax_wendroff_3d(concentration, c_vol, flow, porosity, wells);

    for(auto k = 0; k < n_z - 1; ++k)
        for(auto j = 0; j < n_y - 1; ++j)
            for(auto i = 0; i < n_x - 1; ++i)
            {
                auto index = i + j * n_x + k * n_x * n_y;
                temp_val[index] = c_vol[index].center;
            }

    data_saver.add_step(0, save_data);

    prof.toc("time cycle");

    std::cout << prof;

    for(auto q = ctx.queue().begin(); q != ctx.queue().end(); ++q)
        q->finish();

    return 0;
}
