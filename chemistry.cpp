#include <vector>
#include "parameters.hpp"
#include "utils.hpp"
#include "chemistry.hpp"

//get flux on edges of the cell
inline double edge_center(std::vector<double> &val, int index, int d1)
{
    return 0.5 * (val[index] + val[index + d1]);
}

//get flux through the face
inline double face_center(std::vector<double> &val, int index, int d1, int d2)
{
    return 0.25 * (val[index] + val[index + d1]  + val[index + d2] + val[index + d1 + d2]);
}

//get flux in the cell center
inline double cell_center(std::vector<double> &val, int index)
{
    return 0.125 * (
					val[index] 
					+ val[index + h_i]  
					+ val[index + h_j] 
					+ val[index + h_k] 
					+ val[index + h_i + h_j] 
					+ val[index + h_i + h_k] 
					+ val[index + h_j + h_k] 
					+ val[index + h_i + h_j + h_k]
				   );                                                            
}

//calculate fluxes on the faces of a cell
int get_conc_face(        
					std::vector<double> &face,         
					std::vector<double> &c,
					int idx, 
					int d1, 
					int d2, 
					std::vector<double> &p,
					std::vector<double> &K,
					std::vector<double> &phi
				 )
{
    double tmp;
    double p1, p2, q1, q2;
    double K_tmp = face_center(K, idx, d1, d2);
    double phi_tmp = face_center(phi, idx, d1, d2);

    //clockwise
    //d2 is "senior" axis and heads upward
    //d1 is "junior" axis and heads right
    //chain of command: x-axis < y-axis < z-axis
    //i,j,k=index - lower left corner

    p2 = edge_center(p, idx + d2, d1);
    p1 = edge_center(p, idx, d1);
    q1 = K_tmp * undim * (p2 - p1) / (h * phi_tmp);

    p2 = edge_center(p, idx + d1, d2);
    p1 = edge_center(p, idx, d2);
    q2 = K_tmp * undim * (p2 - p1) / (h * phi_tmp);

    //up - 12 o'clock
    face[0] = edge_center(c, idx + d2, d1);

    //right - 3 o'clock
    face[1] = edge_center(c, idx + d1, d2);

    //down - 6 o'clock
    face[2] = edge_center(c, idx, d1);

    //left - 9 o'clock
    face[3] = edge_center(c, idx, d2);

    //center
    tmp = q2 * ( face[1] - face[3] ) + q1 * ( face[0] - face[2] );
    face[4] = face_center(c, idx, d1, d2) + h_t * tmp / (4. * h);

    return 0;
}

//calculate volume concentration
Cell get_conc_cell(
					std::vector<double> &c, 
					int idx, 
					std::vector<double> &p,
					std::vector<double> &K,
					std::vector<double> &phi
				  )			
{
    Cell tmp;

    //i
    get_conc_face(tmp.x_left, c, idx, h_j, h_k, p, K, phi);

    //i+1
    get_conc_face(tmp.x_right, c, idx + h_i, h_j, h_k, p, K, phi);

    //j
    get_conc_face(tmp.y_left, c, idx, h_i, h_k, p, K, phi);

    //j+1
    get_conc_face(tmp.y_right, c, idx + h_j, h_i, h_k, p, K, phi);

    //k
    get_conc_face(tmp.z_left, c, idx, h_i, h_j, p, K, phi);

    //k+1
    get_conc_face(tmp.z_right, c, idx + h_k, h_i, h_j, p, K, phi);               

    double tmp_x, tmp_y, tmp_z; 
    double p1, p2;
    double K_tmp = cell_center(K, idx); 
    double phi_tmp = cell_center(phi, idx); 

    p2 = face_center(p, idx + h_i, h_j, h_k);
    p1 = face_center(p, idx, h_j, h_k);
    tmp_x = undim * ( p2 - p1 ) * ( tmp.x_right[4] - tmp.x_left[4] ) / h;
 
    p2 = face_center(p, idx + h_j, h_i, h_k);
    p1 = face_center(p, idx, h_i, h_k);
    tmp_y = undim * ( p2 - p1 ) * ( tmp.y_right[4] - tmp.y_left[4] ) / h;

    p2 = face_center(p, idx + h_k, h_i, h_j);
    p1 = face_center(p, idx, h_i, h_j);
    tmp_z = undim * ( p2 - p1 ) * ( tmp.z_right[4] - tmp.z_left[4] ) / h;
   
    tmp.center = cell_center(c, idx) + K_tmp * h_t * (tmp_x + tmp_y + tmp_z) / (2. * h * phi_tmp);

    return tmp;
}

//fill vector of volume concentration
std::vector<Cell> &c_vol get_c_vol(
									std::vector<double> &c, 
									std::vector<double> &p,
									std::vector<double> &K,
									std::vector<double> &phi
								  )
{
	std::vector<Cell> tmp(n);
    for(auto idx = 0; idx < n; ++idx)
		if(!is_border(idx)
			c_vol[idx] = get_conc_cell(c, idx, p, K, phi);

    return tmp;
}

//calculate flux of volume concentration
std::vector<double> get_c_flux(
								std::vector<cell> &c_vol, 
								std::vector<double> &p, 
								std::vector<double> &K, 
								std::vector<double> &phi, 
								int d1, 
								int d2,
								int d3
							  )
{
	std::vector<double> flux(n, 0);

    for(auto idx = 0; idx < n; ++idx)
		if(!is_border(idx)
                flux[idx] = 0.25 *  (
									c_vol[idx].center 
									+ c_vol[idx - d1].center 
									+ c_vol[idx - d2].center
									+ c_vol[idx - d1 - d2].center
									);
}

//solves reactive transport with Corrected Friedrichs scheme for 3D from Liska et al.
void corrected_friedrichs(Unknowns &state)
{
    double tmp;
	//fluxes
    std::vector<double> F(n, 0.);
    std::vector<double> G(n, 0.);
    std::vector<double> V(n, 0.);
	
    double tmp_x, tmp_y, tmp_z;
    double p_tmp;

    state.c_vol = get_c_vol(state.concentration, state.pressure, state.permeability, state.porosity);
    
    F = get_c_flux(state.c_vol, state.pressure, state.permeability, state.porosity, h_j, h_k, h_i);
    G = get_c_flux(state.c_vol, state.pressure, state.permeability, state.porosity, h_i, h_k, h_j);
    V = get_c_flux(state.c_vol, state.pressure, state.permeability, state.porosity, h_i, h_j, h_k);

    for(auto idx = 0; idx < n; ++idx)
		if(!is_border(idx)
        {
			p_tmp = ( p[idx + h_i] - p[idx - h_i] ) / ( 2. * h );
			tmp_x = ( F[idx + h_i] - F[idx - h_i] ) * undim * p_tmp * state.permeability[idx] / state.porosity[idx];

			p_tmp = ( p[idx + h_j] - p[idx - h_j] ) / ( 2. * h );
			tmp_y = ( G[idx] - G[idx - h_j] ) * undim * p_tmp * state.permeability[idx] / state.porosity[idx];

			p_tmp = ( p[idx + h_k] - p[idx - h_k] ) / ( 2. * h );
			tmp_z = ( V[idx] - V[idx - h_k] ) * undim * p_tmp * state.permeability[idx] / state.porosity[idx];

			tmp = ( tmp_x + tmp_y + tmp_z ) * h_t / h;
			
			if(!is_well(idx, wells))
				state.concentration[idx] += tmp;
			else
				state.concentration[idx] = 0.;
        }
}

//update vectors of source, permeability, porosity
void update_s_K_phi(Unknowns &state)
{
#pragma omp parallel for schedule(static)
    for(int index = 0; index < n; ++index)
    {
        if( !(is_well(index, state.wells)) )
        {
            double tearoff = (state.flow[index] < q_0 ? 0. : alfa);

            state.source[index] = beta * state.concentration[index] * ( 1. + b * (phi[idx] - fi_0) )
                - tearoff * state.flow[index] / ( state.porosity[index] * ro_s ) 
                + tearoff * q_0 / ro_s;

            state.porosity[index] += 
                h_t * ( 
                        ( 1. - state.porosity[index] ) * state.dil_dt[index] 
                        - time_un * state.source[index] 
                      ); 

            if(state.porosity[index] > fi_max)
                state.porosity[index] = fi_max;
			
			state.permeability[idx] = koz_car_coef * state.porosity[idx] * state.porosity[idx] * state.porosity[idx] 
				/ ( ( 1. - state.porosity[idx] ) * ( 1. - state.porosity[idx] ) );
        }
    }
}

void solve_chemistry(Unknowns &state)
{
	corrected_friedrichs(state);
	update_s_K_phi(state);
}