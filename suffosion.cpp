#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "boost/multi_array.hpp"
#include <vexcl/vexcl.hpp>

#ifdef nDEBUG
#  undef nDEBUG
#endif
#include <cassert>

/*AMGCL include section*/
#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/io/mm.hpp>
/*end of section*/

#include "pr_util.h"
#include "omp.h"

namespace amgcl
{
    profiler<> prof("project suffosion");
}
using amgcl::prof;

int main(int argc, char *argv[])
{
    // Initialize VexCL context.
    vex::Context ctx( vex::Filter::Env && vex::Filter::DoublePrecision );

    if (!ctx)
    {
        std::cerr << "no GPUs" << std::endl;
        return 1;
    }

    std::cout << ctx << std::endl;

    //Vectors of variables
    std::vector<double> pressure(n, p_top);
    vec_3d velocity(n); //velocities, x,y,z components
    std::vector<double> source(n, 0.);//source of solid phase
    std::vector<double> concentration(n, c_0); 
    std::vector<double> tmp_conc(n, c_0); 
    std::vector<double> permeability(n, 0.);
    std::vector<double> matrix(n, 1.); //matrix of permeability
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
    vex::vector<double> rhs_dev_u(ctx.queue(), rhs_u); // rhs on device
    vex::vector<double> x_ux(ctx.queue(), u_x); // unknown vector on device
    vex::vector<double> x_uy(ctx.queue(), u_y); // unknown vector on device
    vex::vector<double> x_uz(ctx.queue(), u_z); // unknown vector on device
    vex::vector<double> x_u(ctx.queue(), disp);

    //other variables
    const auto duration = ( (argc > 1) ? std::stoul( argv[1] ) : 1 ) * n_t; 
    int writer_step = 0;

    add_well(n_x/2, n_y/2, wells);

    std::sort(wells.begin(), wells.end());

    for(auto i = 0; i < n; ++i)
    {
        double s = 6 * ( 1 - porosity[i] ) / d;
        permeability[i] = pow(porosity[i], 3) / ( T * T * s * s * eta);
    }

    // Define the AMG type:
    typedef amgcl::make_solver<
        amgcl::amg<
        amgcl::backend::vexcl<double>,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
        //amgcl::coarsening::ruge_stuben,
        //amgcl::relaxation::damped_jacobi
            >,
        amgcl::solver::gmres<
            amgcl::backend::vexcl<double>
            >
            > Solver;

    col_pr.reserve(7 * n);
    val_pr.reserve(7 * n);
    ptr_pr.reserve(n);

    prof.tic("build press matrix");
    build_press_mat(col_pr, val_pr, ptr_pr, rhs_pr, permeability, wells);
    Solver solve( boost::tie(n, ptr_pr, col_pr, val_pr) );
    prof.toc("build press matrix");

    prof.tic("solve press matrix");
    vex::copy(rhs_pr, rhs_dev_pr);
    solve(rhs_dev_pr, x_pr);
    vex::copy(x_pr, pressure);   
    prof.toc("solve press matrix");
    
    int n_u = 3 * n;
    col_u.reserve(5 * n_u);
    val_u.reserve(5 * n_u);
    ptr_u.reserve(n_u);
    int    iters;
    double error;

    prof.tic("build disp matrix");
    build_disp_mat(col_u, val_u, ptr_u, wells);
    fill_disp_rhs(pressure, rhs_u, wells);
    Solver solve_disp( boost::tie(n_u, ptr_u, col_u, val_u) );
    prof.toc("build disp matrix");

    prof.tic("solve disp");
    vex::copy(rhs_u, rhs_dev_u);
    boost::tie(iters, error) = solve_disp(rhs_dev_u, x_u);
    vex::copy(x_u, disp);
    prof.toc("solve disp");

    std::cout << "iters = " << iters << std::endl;
    std::cout << "error = " << error << std::endl;
        
    dil_calc(disp, dilatation, dil_dt);
    std::fill(dil_dt.begin(), dil_dt.end(), 0.);

    prof.tic("time cycle");
    for(auto t = 0; t < 0; ++t)
    {
        //pressure
        build_press_mat(col_pr, val_pr, ptr_pr, rhs_pr,permeability, wells);
        vex::copy(rhs_pr, rhs_dev_pr);
        Solver solve(boost::tie(n, ptr_pr, col_pr, val_pr));
        solve(rhs_dev_pr, x_pr);
        vex::copy(x_pr, pressure);

        //velocity
        vel_calc(pressure, velocity, permeability, wells);

        //concentration
        tmp_conc = concentration;
        conc_calc(concentration, tmp_conc, porosity, source, velocity, wells);

        //displacements
        fill_disp_rhs(pressure, rhs_u, wells);
        vex::copy(rhs_u, rhs_dev_u);
        vex::copy(x_u, disp);
        
        //dilatation
        dil_calc(disp, dilatation, dil_dt);

        //porosity&source
        por_calc(concentration, porosity, source, velocity, dil_dt);

        //permeability
        per_calc(porosity, permeability);

        //write results
        /*       writer_step++;
                 if( writer_step == 100 )
                 {
                 std::string filename = "./data/pressure" + std::to_string(t) + ".vtk";
                 wrt_vtk(pressure, filename);

                 filename = "./data/qx" + std::to_string(t) + ".vtk";
                 wrt_vtk(velocity.x, filename);

                 filename = "./data/qy" + std::to_string(t) + ".vtk";
                 wrt_vtk(velocity.y, filename);

                 filename = "./data/qz" + std::to_string(t) + ".vtk";
                 wrt_vtk(velocity.z, filename);

                 filename = "./data/dil" + std::to_string(t) + ".vtk";
                 wrt_vtk(dilatation, filename);

                 filename = "./data/porosity" + std::to_string(t) + ".vtk";
                 wrt_vtk(porosity, filename);

                 filename = "./data/dphi_dt" + std::to_string(t) + ".vtk";
                 wrt_vtk(por_dt, filename);

                 filename = "./data/dil_dt" + std::to_string(t) + ".vtk";
                 wrt_vtk(dil_dt, filename);

                 filename = "./data/permeability" + std::to_string(t) + ".vtk";
                 wrt_vtk(permeability, filename);

                 filename = "./data/concentration" + std::to_string(t) + ".vtk";
                 wrt_vtk(concentration, filename);

                 writer_step = 0;

                 }*/
    }
    
    for(auto index = 0; index < n; ++index)
    {
        u_x[index] = disp[3*index + 0];
        u_y[index] = disp[3*index + 1];
        u_z[index] = disp[3*index + 2];
    }

    prof.toc("time cycle");

    std::cout << prof;

    for(auto q = ctx.queue().begin(); q != ctx.queue().end(); ++q)
        q->finish();

    std::string filename = "./data/pressure_fin.vtk";
    wrt_vtk(pressure, filename);

    filename = "./data/qx_fin.vtk";
    wrt_vtk(velocity.x, filename);

    filename = "./data/qy_fin.vtk";
    wrt_vtk(velocity.y, filename);

    filename = "./data/qz_fin.vtk";
    wrt_vtk(velocity.z, filename);

    filename = "./data/ux_fin.vtk";
    wrt_vtk(u_x, filename);

    filename = "./data/uy_fin.vtk";
    wrt_vtk(u_y, filename);

    filename = "./data/uz_fin.vtk";
    wrt_vtk(u_z, filename);

    filename = "./data/dil_fin.vtk";
    wrt_vtk(dilatation, filename);

    filename = "./data/porosity_fin.vtk";
    wrt_vtk(porosity, filename);

    filename = "./data/por_dt_fin.vtk";
    wrt_vtk(por_dt, filename);

    filename = "./data/dil_dt_fin.vtk";
    wrt_vtk(dil_dt, filename);

    filename = "./data/permeability_fin.vtk";
    wrt_vtk(permeability, filename);

    filename = "./data/concentration_fin.vtk";
    wrt_vtk(concentration, filename);

    filename = "./data/source_fin.vtk";
    wrt_vtk(source, filename);

    return 0;
}
