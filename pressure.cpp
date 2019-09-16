#include <cmath>
#include <tuple>
#include "pressure.hpp"
#include "parameters.hpp"

/*AMGCL include section*/
#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/gmres.hpp>
/*end of section*/

//calculate abs value of 3d flow
std::vector<double> get_flow(std::vector<double> &p, std::vector<double> &K)
{
    double qx = 0., qy = 0., qz = 0.;
    const double loc_undim = lame_2 / length ;
	std::vector<double> q(n);

    boost::array<char, 3> border_flag = {0, 0, 0};

    for(auto idx = 0; idx < n; ++idx)
    {
        border_flag = check_border(idx);

        //x-axis flow
        if(border_flag[0] != 0)
            if(border_flag[0] == 'w')
                qx = loc_undim * K[idx] * (p[idx + h_i] - p[idx]) / h;
            else       
                qx = loc_undim * K[idx] * (p[idx] - p[idx - h_i]) / h;
        else
            qx = loc_undim * K[idx] * (p[idx + h_i] - p[idx - h_i]) / (2. * h);

        //y-axis flow
        if(border_flag[1] != 0)
            qx = qy = qz = 0.;
        else
            qy = loc_undim * K[idx] * (p[idx + h_j] - p[idx - h_j]) / (2. * h);

        //z-axis flow
        if(border_flag[2] != 0)
            qx = qy = qz = 0.;
        else
            qz = loc_undim * K[idx] * (p[idx + h_k] - p[idx - h_k]) / (2. * h);

        q[idx] = sqrt( (qx * qx) + (qy * qy) + (qz * qz) );
    }
    return q;    
}

//solve pressure equation
void solve_pressure(Unknowns &state, vex::Context &ctx, amgcl::profiler &prof)
{
    //set backend
    typedef amgcl::backend::vexcl<double>   SBackend;
    
    // Define the AMG type:
    typedef amgcl::make_solver<
        amgcl::amg<
        SBackend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
            >,
        amgcl::solver::gmres< SBackend >
            > PSolver;
	
	//Vectors for calculations
    std::vector<double> val; //values of nonzero entries
    std::vector<int>    col; //column numbers of nonzero entries
    std::vector<int>    ptr; //points to the start of each row in col
    std::vector<double> rhs(n, 0); //right-hand side 
    vex::vector<double> rhs_dev(ctx.queue(), rhs_pr); //rhs on device
    vex::vector<double> x(ctx.queue(), state.pressure); //unknown vector on device
	
    col.reserve(7 * n);
    val.reserve(7 * n);
    ptr.reserve(n);
	
	int    iters;
    double error;

    prof.tic("build pressure matrix");
    assembly_matrix(col, val, ptr, state);
    PSolver solve( std::tie(n, ptr, col, val) );
    prof.toc("build pressure matrix");

    prof.tic("solve pressure");
    vex::copy(rhs, rhs_dev);
    std::tie(iters, error) = solve(rhs_dev, x);
    vex::copy(x, state.pressure);   
    prof.toc("solve pressure");

    std::cout << "iters P = " << iters << std::endl;
    std::cout << "error P = " << error << std::endl;
	
	state.flow = get_flow(state.pressure, state.permeability);
}

//consider log-distribution of pressure around well
double get_nu(double i, double j, std::set<std::pair<int, int>> wells)
{
    double sum = 0.;

    for(int index = 0; index < wells.size(); index += n_z)
    {
        int z = index / (n_x * n_y);
        int tmp_loc = index - z * n_x * n_y;
        int y = tmp_loc / n_x;
        int x = tmp_loc - y * n_x;

        double r_2 = pow(double(x) - i, 2.) + pow(double(y) - j, 2.);

        if( fabs(r_2 - 0.25) < 0.0001 )
            sum += 3.14 / ( 2.0 * log(h/r_c) );
        else if( fabs(r_2 - 1.25) < 0.0001 )
            sum += 2.0 * atan(0.5) / log(2.0);
        else if( fabs(r_2 - 2.25) < 0.0001 )
            sum += 2.0 * atan(1.0/3.0) / log(2.0);
    }

    if(fabs(sum) < 1e-7)
        sum = 1.;

    return sum;
}

//calculate mobilities
double get_coefficient(int index_1, int index_2, std::vector<double> &K, std::set<std::pair<int, int>> &wells)
{
    int k_1 = index_1 / (n_x * n_y);
    int tmp_loc = index_1 - k_1 * n_x * n_y;
    int j_1 = tmp_loc / n_x;
    int i_1 = tmp_loc - j_1 * n_x;

    int k_2 = index_2 / (n_x * n_y);
    tmp_loc = index_2 - k_2 * n_x * n_y;
    int j_2 = tmp_loc / n_x;
    int i_2 = tmp_loc - j_2 * n_x;

    double perm = ( K[index_1] + K[index_2] ) / 2.;

    double tmp = perm * time_un * lame_2 / (length * length);
 
    int dx = i_2 - i_1;
    int dy = j_2 - j_1;

    if(dx)
        tmp *= get_nu(double(i_1) + dx * 0.5, double(j_1), wells);
     
    if(dy)
        tmp *= get_nu(double(i_1), double(j_1) + dy * 0.5, wells);
	
    return -tmp;
}

//build matrix for pressure equation
void assembly_matrix(
						std::vector<int> &col, 
						std::vector<double> &val, 
						std::vector<int> &ptr, 
						std::vector<double> &rhs, 
						Unknowns &state
					)
{    
    std::fill(rhs.begin(), rhs.end(), 0.);

    col.clear();  ptr.clear(); val.clear();  
    
    ptr.push_back(0);

    boost::array<char, 3> border_flag = {0, 0, 0};

    for(int index = 0; index < n; ++index)
    {
        rhs[index] = h * h * (state.source[index] - state.dil_dt[index]);

		//7-stencil coefficients
        double k_f, k_b, j_f, j_b, i_f, i_b, cntr;

        border_flag = check_border(index);

		//check if boundary
        if(border_flag[0] != 0)
        {
            if(border_flag[0] == 'w')
            {
                rhs[index] = p_bh;
                cntr = 1.;
                k_f = k_b = j_f = j_b = i_b = i_f = 0.;
            }
            else
            {
                rhs[index] = p_top;
                cntr = 1.;
                k_f = k_b = j_f = j_b = i_b = i_f = 0.;
            }
        }
        else
        {
            i_f = get_coefficient(index, index + h_i, state.permeability, state.wells);
            i_b = get_coefficient(index, index - h_i, state.permeability, state.wells);
        }               

        if(border_flag[1] != 0)
        {
            if(border_flag[1] == 's')
                j_f = 2. * get_coefficient(index, index + h_j, state.permeability, state.wells);
            else
                j_b = 2. * get_coefficient(index, index - h_j, state.permeability, state.wells);        
        }
        else
        {
            j_f = get_coefficient(index, index + h_j, state.permeability, state.wells);
            j_b = get_coefficient(index, index - h_j, state.permeability, state.wells);
        }

        if(border_flag[2] != 0)
        {
            if(border_flag[2] == 'b')
                k_f = 2. * get_coefficient(index, index + h_k, state.permeability, state.wells);
            else
                k_b = 2. * get_coefficient(index, index - h_k, state.permeability, state.wells);
        }
        else
        {
            k_f = get_coefficient(index, index + h_k, state.permeability, state.wells);
            k_b = get_coefficient(index, index - h_k, state.permeability, state.wells);
        }

        cntr += -(k_f + k_b + j_b + j_f + i_b + i_f);
		
		//check if well
        if(is_well(index, wells))
        {
            rhs[index] = p_bh;
            cntr = 1.;
            k_f = k_b = j_f = j_b = i_b = i_f = 0.;
        }

        if(is_well(index + h_j, state.wells))
        {
            rhs[index] -= j_f * p_bh;
            j_f = 0.;
        }

        if(is_well(index - h_j, state.wells))
        {
            rhs[index] -= j_b * p_bh;
            j_b = 0.;
        }

        if(is_well(index + h_i, state.wells))
        {
            rhs[index] -= i_f * p_bh;
            i_f = 0.;
        }

        if(is_well(index - h_i, state.wells))
        {
            rhs[index] -= i_b * p_bh;
            i_b = 0.;
        }

        place_val(k_b, index - n_x * n_y, col, val);
        place_val(j_b, index - n_x, col, val);
        place_val(i_b, index - 1, col, val);
        place_val(cntr, index, col, val);
        place_val(i_f, index + 1, col, val);
        place_val(j_f, index + n_x, col, val);
        place_val(k_f, index + n_x * n_y, col, val);

        ptr.push_back(col.size());
    }
}