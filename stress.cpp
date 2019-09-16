#include <tuple>
#include "parameters.hpp"
#include "utils.hpp"
#include "stress.hpp"

/*AMGCL-VexCL include section*/
#include <amgcl/amg.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/backend/vexcl_static_matrix.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/solver/bicgstab.hpp>
/*end of section*/

//update right-hand for Lame system
std::vector<double> fill_rhs(std::vector<double> &pressure, std::set<std::pair<int, int>> &wells)
{
    double tmp = 0.;
	std::vector<double> rhs(3 * n, 0);
	
    for(int index = 0; index < n; ++index)
    {
        if( !(is_well(index, wells)) && !(is_border(index)) )
        {
            if(pressure[index] < pressure[index - h_i])
                tmp = (pressure[index] - pressure[index - h_i]);
            else
                tmp = (pressure[index + h_i] - pressure[index]);
            rhs[3*index + 0] = -h * tmp;

            if(pressure[index] < pressure[index - h_j])
                tmp = (pressure[index] - pressure[index - h_j]);
            else
                tmp = (pressure[index + h_j] - pressure[index]);
            rhs[3*index + 1] = -h * tmp;

            if(pressure[index] < pressure[index - h_k])
                tmp = (pressure[index] - pressure[index - h_k]);
            else
                tmp = (pressure[index + h_k] - pressure[index]);               
            rhs[3*index + 2] = -h * tmp;

        }
        else
            rhs[index] = 0.;
    }

    return rhs;
}

//solve Lame system
void solve_stress(Unknowns &state, vex::Context &ctx, amgcl::profiler &prof)
{
    typedef amgcl::static_matrix<double, 3, 3> mat_type;
    typedef amgcl::static_matrix<double, 3, 1> vec_type;

    vex::scoped_program_header header(ctx, amgcl::backend::vexcl_static_matrix_declaration<double,3>());
	
    //Vectors for calculations
    std::vector<double> val; //values of nonzero entries
    std::vector<int>    col; //column numbers of nonzero entries
    std::vector<int>    ptr; //points to the start of each row in col_pr 
    std::vector<double> rhs(3 * n, 0); // right-hand side 
    std::vector<double> disp(3 * n, 0); //block displacements
    vex::vector<vec_type> rhs_dev(ctx.queue(), n); // rhs on device
    vex::vector<vec_type> x(ctx.queue(), n);
    amgcl::backend::clear(x);	
	   
	// Define the AMG type:
    typedef amgcl::backend::vexcl<mat_type> BBackend;

    typedef amgcl::make_solver<
        amgcl::amg<
        BBackend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::damped_jacobi
            >,
        amgcl::solver::bicgstab< BBackend >
            > USolver;
	
	int n_u = 3 * n;
    col.reserve(5 * n_u);
    val.reserve(5 * n_u);
    ptr.reserve(n_u);
	
	prof.tic("build Lame matrix");
    assembly_matrix(col, val, ptr, rhs, state.pressure, state.wells);
    USolver solve( 
					amgcl::adapter::block_matrix<mat_type>
					(
						boost::tie(n_u, ptr_u, col_u, val_u) 
					) 
				 );
    prof.toc("build Lame matrix");

    prof.tic("solve disp");
    vec_type const * fptr = reinterpret_cast<vec_type const *>(&rhs[0]);
    vec_type       * xptr = reinterpret_cast<vec_type       *>(&disp[0]);
    vex::copy(fptr, fptr + n, rhs_dev.begin());
    std::tie(iters, error) = solve_disp(rhs_dev, x);
    vex::copy(x.begin(), x.end(), xptr);
    prof.toc("solve disp");

    std::cout << "iters U = " << iters << std::endl;
    std::cout << "error U = " << error << std::endl;
	
    for(auto index = 0; index < n; ++index)
    {
        state.u_x[index] = disp[3 * index + 0];
        state.u_y[index] = disp[3 * index + 1];
        state.u_z[index] = disp[3 * index + 2];
    }	

    get_dilatation(state);
}

//build matrix for Lame system
void assembly_matrix(
						std::vector<int> &col, 
						std::vector<double> &val, 
						std::vector<int> &ptr,
						std::vector<double> &rhs,						
						std::vector<double> &pressure, 
						std::set<std::pair<int, int>> wells
					)
{
    double c_1 = 1. + lame_1/lame_2;
    double c_2= 1. + c_1;          

    ptr.push_back(0);
    int point;

    for(int index = 0; index < n; ++index)
    {
		//fill rhs
		rhs = fill_rhs(pressure, wells);
		
		//check border
        if (is_border(index) || is_well(index, wells)) 
        {
            place_val(1., X(index), col, val);
            ptr.push_back(col.size());

            place_val(1., Y(index), col, val);
            ptr.push_back(col.size());

            place_val(1., Z(index), col, val);
            ptr.push_back(col.size());
        }
        else
        {
            //u_x line

            //i,j,k-1
            point = index - h_k;
            //u_x
            place_val(-1., X(point), col, val);
            //u_z
            place_val(-c_1, Z(point), col, val);

            //i+1,j,k-1
            point = index - h_k + h_i;
            if(!is_well(point, wells))
                //u_z
                place_val(c_1, Z(point), col, val);

            //i,j-1,k
            point = index - h_j;
            if(!is_well(point, wells))
            {
                //u_x
                place_val(-1., X(point), col, val);
                //u_y
                place_val(-c_1, Y(point), col, val);
            }

            //i+1,j-1,k
            point = index - h_j + h_i;
            if(!is_well(point, wells))
                //u_y
                place_val(c_1, Y(point), col, val);

            //i-1,j,k
            point = index - h_i;
            if(!is_well(point, wells))
                //u_x
                place_val(-c_2, X(point), col, val);

            //i,j,k
            point = index;
            //u_x            
            place_val(2. * c_2 + 4., X(point), col, val);
            //u_y
            place_val(c_1, Y(point), col, val);
            //u_z
            place_val(c_1, Z(point), col, val);

            //i+1,j,k
            point = index + h_i;
            if(!is_well(point, wells))
            {
                //u_x            
                place_val(-c_2, X(point), col, val);
                //u_y
                place_val(-c_1, Y(point), col, val);
                //u_z
                place_val(-c_1, Z(point), col, val);
            }

            //i,j+1,k
            point = index + h_j;
            if(!is_well(point, wells))
                //u_x
                place_val(-1., X(point), col, val);

            //i,j,k+1
            point = index + h_k;
            //u_x
            place_val(-1., X(point), col, val);

            ptr.push_back(col.size());


            //u_y line

            //i,j,k-1
            point = index - h_k;
            //u_y
            place_val(-1., Y(point), col, val);

            //i,j-1,k
            point = index - h_j;
            if(!is_well(point, wells))
            {
                //u_x
                place_val(-c_1, X(point), col, val);         
                //u_y
                place_val(-c_2, Y(point), col, val);         
                //u_z
                place_val(-c_1, Z(point), col, val);         
            }

            //i+1,j-1,k
            point = index - h_j + h_i;
            if(!is_well(point, wells))
                //u_x
                place_val(c_1, X(point), col, val);         

            //i-1,j,k
            point = index - h_i;
            if(!is_well(point, wells))
                //u_y
                place_val(-1., Y(point), col, val);         

            //i,j,k
            point = index;
            //u_x            
            place_val(c_1, X(point), col, val);
            //u_y
            place_val(2. * c_2 + 4., Y(point), col, val);
            //u_z
            place_val(c_1, Z(point), col, val);           

            //i+1,j,k
            point = index + h_i;
            if(!is_well(point, wells))
            {
                //u_x
                place_val(-c_1, X(point), col, val);
                //u_y
                place_val(-1., Y(point), col, val);         
            }

            //i,j+1,k
            point = index + h_j;
            //u_y
            if(!is_well(point, wells))
                //u_y
                place_val(-c_2, Y(point), col, val);

            //i,j-1,k+1
            point = index - h_j + h_k;
            if(!is_well(point, wells))
                //u_z
                place_val(c_1, Z(point), col, val);

            //i,j,k+1
            point = index + h_k;
            //u_y
            place_val(-1., Y(point), col, val);
            //u_z
            place_val(-c_1, Z(point), col, val);

            ptr.push_back(col.size());


            //u_z line

            //i,j,k-1
            point = index - h_k;
            //u_x
            place_val(-c_1, X(point), col, val);
            //u_z
            place_val(-c_2, Z(point), col, val);

            //i+1,j,k-1
            point = index - h_k + h_i;
            if(!is_well(point, wells))
                //u_x
                place_val(c_1, X(point), col, val);

            //i,j-1,k
            point = index - h_j;
            if(!is_well(point, wells))
            {
                //u_y
                place_val(-c_1, Y(point), col, val);
                //u_z
                place_val(-1., Z(point), col, val);
            }

            //i-1,j,k
            point = index - h_i;
            if(!is_well(point, wells))
                //u_z
                place_val(-1., Z(point), col, val);

            //i,j,k
            point = index;
            //u_x            
            place_val(c_1, X(point), col, val);
            //u_y
            place_val(c_1, Y(point), col, val);
            //u_z
            place_val(2. * c_2 + 4., Z(point), col, val);

            //i+1,j,k
            point = index + h_i;
            if(!is_well(point, wells))
            {
                //u_x
                place_val(-c_1, X(point), col, val);
                //u_z
                place_val(-1., Z(point), col, val);
            }

            //i,j+1,k
            point = index + h_j;
            //u_z
            if(!is_well(point, wells))
                //u_z
                place_val(-1., Z(point), col, val);

            //i,j-1,k+1
            point = index - h_j + h_k;
            if(!is_well(point, wells))
                //u_y
                place_val(c_1, Y(point), col, val);

            //i,j,k+1
            point = index + h_k;
            //u_y
            place_val(-c_1, Y(point), col, val);         
            //u_z
            place_val(-c_2, Z(point), col, val);         

            ptr.push_back(col.size());
        }
    }
}

//calculate trace of strain tensor
void get_dilatation(Unknowns &state)
{
    double tmp = 0.;

    for(int index = 0; index < n; ++index)
    {
        tmp = 0.;

		//check boundary
        if( !(is_border(index)) && !(is_well(index, wells)) )
        {
            if(!(is_well(index + h_i, wells)) && !(is_well(index - h_i, wells)))
                tmp += u_x(index + h_i) - u_x(index - h_i);
            
            if(!(is_well(index + h_j, wells)) && !(is_well(index - h_j, wells)))
                tmp += u_y(index + h_j)] - u_y(index - h_j);
    
            tmp += u_z(index + h_k) - u_z(index - h_k);

            tmp /= (2. * h);

            state.dil_dt[index] = ( tmp - state.dilatation[index] ) / (h_t * time_un);

            state.dilatation[index] = tmp;
        }
        else
        {
            state.dilatation[index] = 0.;
            state.dil_dt[index] = 0.;
        }
    }
}
