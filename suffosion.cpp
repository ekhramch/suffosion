#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
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
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/profiler.hpp>
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
    std::vector<double> pressure(n, 0.);
    vec_3d velocity(n); //velocities, x,y,z components
    std::vector<double> source(n, 0.);//source of solid phase
    std::vector<double> concentration(n, c_0); 
    std::vector<double> permeability(n, 0.);
    std::vector<double> matrix(n, 1.); //matrix of permeability
    std::vector<double> u_x(n, 0.), u_y(n, 0.), u_z(n, 0.); //displacements
    std::vector<double> dilatation(n, 0.); 
    std::vector<double> dil_dt(n, 0.); //dilatation time derivative
    std::vector<double> porosity(n, fi_0);
    std::vector<double> por_dt(n, 0.); //porosity time derivative
    std::vector<point> boreholes; //coordinates of boreholes

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
    std::vector<double> rhs_u(n, 0.); // right-hand side 
    vex::vector<double> rhs_dev_u(ctx.queue(), rhs_u); // rhs on device
    vex::vector<double> x_ux(ctx.queue(), u_x); // unknown vector on device
    vex::vector<double> x_uy(ctx.queue(), u_y); // unknown vector on device
    vex::vector<double> x_uz(ctx.queue(), u_z); // unknown vector on device

    //other variables
    const auto time = ( (argc > 1) ? std::stoul( argv[1] ) : 1 ) * n_t; 
    int writer_step = 0;

    for(auto i = 0; i < n; ++i)
    {
        double s = 6 * ( 1 - porosity[i] ) / d;
        permeability[i] = pow(porosity[i], 3) / ( T * T * s * s * eta);
    }

    // Define the AMG type:
    typedef amgcl::make_solver<
        amgcl::amg<
        amgcl::backend::vexcl<double>,
        amgcl::coarsening::ruge_stuben,
        amgcl::relaxation::damped_jacobi
            >,
        amgcl::solver::gmres<
            amgcl::backend::vexcl<double>
            >
            > Solver;

    prof.tic("build press matrix");
    build_press_mat(col_pr, val_pr, ptr_pr, rhs_pr, permeability, boreholes);
    Solver solve( boost::tie(n, ptr_pr, col_pr, val_pr) );
    prof.toc("build press matrix");

    prof.tic("solve press matrix");
    vex::copy(rhs_pr, rhs_dev_pr);
    solve(rhs_dev_pr, x_pr);
    vex::copy(x_pr, pressure);   
    get_q(pressure, velocity, permeability, boreholes);
    prof.toc("solve press matrix");

    build_u_mat(col_u, val_u, ptr_u, 'x');
    Solver solve_ux( boost::tie(n, ptr_u, col_u, val_u) );
    col_u.clear();  ptr_u.clear(); val_u.clear();

    build_u_mat(col_u, val_u, ptr_u, 'y');
    Solver solve_uy( boost::tie(n, ptr_u, col_u, val_u) );
    col_u.clear();  ptr_u.clear(); val_u.clear();

    build_u_mat(col_u, val_u, ptr_u, 'z');
    Solver solve_uz( boost::tie(n, ptr_u, col_u, val_u) );
    col_u.clear();  ptr_u.clear(); val_u.clear();

    prof.tic("time cycle");
    for(auto t = 0; t < 0; ++t)
    {
        //pressure
/*        std::fill(rhs_pr.begin(), rhs_pr.end(), 0.);
        vex::copy(rhs_pr, rhs_dev_pr);

        build_press_mat(col_pr, val_pr, ptr_pr, rhs_pr,permeability, boreholes);
        Solver solve( boost::tie(n, ptr_pr, col_pr, val_pr) );
        
        vex::copy(rhs_pr, rhs_dev_pr);
        solve(rhs_dev_pr, x_pr);
        
        vex::copy(x_pr, pressure);
        get_q(pressure, velocity, permeability, boreholes);

        //concentration
        conc_calc(concentration, porosity, source, velocity, boreholes, time);

        //u_x
        for(auto index = 0; index < n; ++index)
        {
            auto k = index / h_k;
            auto tmp_loc = index - k * h_k;
            auto j = tmp_loc / h_j;
            auto i = tmp_loc - j * h_j;
            bool is_border =  
                (k == 0) || (k == n_z - 1) || (j == 0) || 
                (j == n_y - 1) || (i == 0) || (i ==  n_x - 1); 

            if(is_border)
                rhs_u[index] = 0.;
            else
                rhs_u[index] = 
                    sec_ord_mx(u_y, index, "xy") +
                    sec_ord_mx(u_z, index, "xz") - 
                    ( pressure[index + h_i] - pressure[index - h_i] ) /
                    ( 2. * lame_2 );
        }
        vex::copy(rhs_u, rhs_dev_u);
        solve_ux(rhs_dev_u, x_ux);
        vex::copy(x_ux, u_x);
        check_disp(u_x);

        //u_y
        for(auto index = 0; index < n; ++index)
        {
            auto k = index / h_k;
            auto tmp_loc = index - k * h_k;
            auto j = tmp_loc / h_j;
            auto i = tmp_loc - j * h_j;
            bool is_border =  
                (k == 0) || (k == n_z - 1) || (j == 0) || 
                (j == n_y - 1) || (i == 0) || (i ==  n_x - 1); 

            if(is_border)
                rhs_u[index] = 0.;
            else
                rhs_u[index] = 
                    sec_ord_mx(u_x, index, "xy") +
                    sec_ord_mx(u_z, index, "yz") - 
                    ( pressure[index + h_j] - pressure[index - h_j] ) /
                    ( 2. * lame_2 );
        }
        vex::copy(rhs_u, rhs_dev_u);
        solve_ux(rhs_dev_u, x_uy);
        vex::copy(x_uy, u_y);
        check_disp(u_y);

        //u_z
        for(auto index = 0; index < n; ++index)
        {
            auto k = index / h_k;
            auto tmp_loc = index - k * h_k;
            auto j = tmp_loc / h_j;
            auto i = tmp_loc - j * h_j;
            bool is_border =  
                (k == 0) || (k == n_z - 1) || (j == 0) || 
                (j == n_y - 1) || (i == 0) || (i ==  n_x - 1); 

            if(is_border)
                rhs_u[index] = 0.;
            else
                rhs_u[index] = 
                    sec_ord_mx(u_x, index, "xz") +
                    sec_ord_mx(u_y, index, "yz") - 
                    ( pressure[index + h_k] - pressure[index - h_k] ) /
                    ( 2. * lame_2 );
        }                   
        vex::copy(rhs_u, rhs_dev_u);
        solve_ux(rhs_dev_u, x_uz);
        vex::copy(x_uz, u_z);
        check_disp(u_z);

        //dilatation
        for (auto k = 1, index = in_idx; k < n_z - 1; ++k)
            for (auto j = 1; j < n_y - 1; ++j)
                for (auto i = 1; i < n_x - 1; ++i, ++ index)
                {
                    double tmp = u_z[index + h_k] - u_z[index - h_k];
                    tmp += u_y[index + h_j] - u_y[index - h_j];
                    tmp += u_x[index + h_i] - u_z[index - h_i];
                    tmp /= (2. * h);

                    if(t > 0)
                        dil_dt[index] = ( tmp - dilatation[index] ) / h_t;

                    dilatation[index] = tmp;
                }

        //porosity&source&permeability
        for(auto index = 0; index < 0; ++index)
        {
            auto vel_mod = sqrt( 
                    pow(velocity.x[index]/q_0, 2.) + 
                    pow(velocity.y[index]/q_0, 2.) +
                    pow(velocity.z[index]/q_0, 2.) );

            double tearoff = (vel_mod < q_0 ? 0. : alfa);

            source[index] = beta * concentration[index] 
                - tearoff * vel_mod / ( porosity[index] * ro_s ) 
                + tearoff * q_0 / ro_s;

            auto tmp = porosity[index];

            porosity[index] += 0.5 * h_t * 
                ( ( 1 - porosity[index] ) * dil_dt[index] - source[index] );

            por_dt[index] = 2. * ( porosity[index] - tmp ) / h_t;

            double s = 6 * ( 1 - porosity[index] ) / d;

            permeability[index] = pow(porosity[index], 3)/(T * T * s * s * eta);
        }
        
        
        //pressure
        for(auto index = 0; index < n; ++index)
            rhs_pr[index] = -h * h *
                ( por_dt[index] + porosity[index] * dil_dt[index] );
        build_press_mat(col_pr, val_pr, ptr_pr, rhs_pr,permeability, boreholes);
        vex::copy(rhs_pr, rhs_dev_pr);
        //solve(rhs_dev_pr, x_pr);
        //vex::copy(x_pr, pressure);
        get_q(pressure, velocity, permeability, boreholes);

        //concentration
        conc_calc(concentration, porosity, source, velocity, boreholes, time);

        //u_x
        for(auto index = 0; index < n; ++index)
        {
            auto k = index / h_k;
            auto tmp_loc = index - k * h_k;
            auto j = tmp_loc / h_j;
            auto i = tmp_loc - j * h_j;
            bool is_border =  
                (k == 0) || (k == n_z - 1) || (j == 0) || 
                (j == n_y - 1) || (i == 0) || (i ==  n_x - 1); 

            if(is_border)
                rhs_u[index] = 0.;
            else
                rhs_u[index] = 
                    sec_ord_mx(u_y, index, "xy") +
                    sec_ord_mx(u_z, index, "xz") - 
                    ( pressure[index + h_i] - pressure[index - h_i] ) /
                    ( 2. * lame_2 );
        }
        vex::copy(rhs_u, rhs_dev_u);
        solve_ux(rhs_dev_u, x_ux);
        vex::copy(x_ux, u_x);
        check_disp(u_x);

        //u_y
        for(auto index = 0; index < n; ++index)
        {
            auto k = index / h_k;
            auto tmp_loc = index - k * h_k;
            auto j = tmp_loc / h_j;
            auto i = tmp_loc - j * h_j;
            bool is_border =  
                (k == 0) || (k == n_z - 1) || (j == 0) || 
                (j == n_y - 1) || (i == 0) || (i ==  n_x - 1); 

            if(is_border)
                rhs_u[index] = 0.;
            else
                rhs_u[index] = 
                    sec_ord_mx(u_x, index, "xy") +
                    sec_ord_mx(u_z, index, "yz") - 
                    ( pressure[index + h_j] - pressure[index - h_j] ) /
                    ( 2. * lame_2 );
        }
        vex::copy(rhs_u, rhs_dev_u);
        solve_ux(rhs_dev_u, x_uy);
        vex::copy(x_uy, u_y);
        check_disp(u_y);

        //u_z
        for(auto index = 0; index < n; ++index)
        {
            auto k = index / h_k;
            auto tmp_loc = index - k * h_k;
            auto j = tmp_loc / h_j;
            auto i = tmp_loc - j * h_j;
            bool is_border =  
                (k == 0) || (k == n_z - 1) || (j == 0) || 
                (j == n_y - 1) || (i == 0) || (i ==  n_x - 1); 

            if(is_border)
                rhs_u[index] = 0.;
            else
                rhs_u[index] = 
                    sec_ord_mx(u_x, index, "xz") +
                    sec_ord_mx(u_y, index, "yz") - 
                    ( pressure[index + h_k] - pressure[index - h_k] ) /
                    ( 2. * lame_2 );
        }                   
        vex::copy(rhs_u, rhs_dev_u);
        solve_ux(rhs_dev_u, x_uz);
        vex::copy(x_uz, u_z);
        check_disp(u_z);

        //dilatation
        for (auto k = 1, index = in_idx; k < n_z - 1; ++k)
            for (auto j = 1; j < n_y - 1; ++j)
                for (auto i = 1; i < n_x - 1; ++i, ++ index)
                {
                    double tmp = u_z[index + h_k] - u_z[index - h_k];
                    tmp += u_y[index + h_j] - u_y[index - h_j];
                    tmp += u_x[index + h_i] - u_z[index - h_i];
                    tmp /= (2. * h);

                    if(t > 0)
                        dil_dt[index] = ( tmp - dilatation[index] ) / h_t;

                    dilatation[index] = tmp;
                }

        //porosity&source&permeability
        for(auto index = 0; index < 0; ++index)
        {
            auto vel_mod = sqrt( 
                    pow(velocity.x[index]/q_0, 2.) + 
                    pow(velocity.y[index]/q_0, 2.) +
                    pow(velocity.z[index]/q_0, 2.) );

            double tearoff = (vel_mod < q_0 ? 0. : alfa);

            source[index] = beta * concentration[index] 
                - tearoff * vel_mod / ( porosity[index] * ro_s ) 
                + tearoff * q_0 / ro_s;

            auto tmp = porosity[index];

            porosity[index] += 0.5 * h_t * 
                ( ( 1 - porosity[index] ) * dil_dt[index] - source[index] );

            por_dt[index] = 2. * ( porosity[index] - tmp ) / h_t;

            double s = 6 * ( 1 - porosity[index] ) / d;

            permeability[index] = pow(porosity[index], 3)/(T * T * s * s * eta);
        }

        //write results
        writer_step++;
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

            filename = "./data/ux" + std::to_string(t) + ".vtk";
            wrt_vtk(u_x, filename);

            filename = "./data/uy" + std::to_string(t) + ".vtk";
            wrt_vtk(u_y, filename);

            filename = "./data/uz" + std::to_string(t) + ".vtk";
            wrt_vtk(u_z, filename);

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
