#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <assert.h>
#include "omp.h"
#include <cmath>
#include <string>
#include "pr_util.h"

#define X(i)(3*(i) + 0)
#define Y(i)(3*(i) + 1)
#define Z(i)(3*(i) + 2)

using namespace std;

int wrt_vtk(vector<double> &arr, const string &filename)
{
    fstream f(filename.c_str(), ios::out);

    f << "# vtk DataFile Version 2.0" << endl;

    f << "data" << endl;

    f << "ASCIIY" << endl;

    f << "DATASET STRUCTURED_POINTS" << endl;

    f << "DIMENSIONS " << n_x << " " << n_y << " " << n_z  << endl;

    f << "ASPECT_RATIO 1.0 1.0 1.0" << endl;

    f << "ORIGIN 0.0 0.0 0.0" << endl;

    f << "POINT_DATA " << n_x * n_y * n_z << endl;

    f << "SCALARS data double" << endl;

    f << "LOOKUP_TABLE default" << endl;

    for(int k = 0; k < n_z; ++k)
    {
	for(int j = 0; j < n_y; ++j)
	    for(int i = 0; i < n_x; ++i)
		f <<  arr[i + j*n_x + k*n_x*n_y] << " ";

	f << endl;
    }

    return 0;
}

double get_nu(double i, double j, vector<int> wells)
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

int conc_calc(vector<double> &conc_1, vector<double> &conc_2, 
        vector<double> porosity, vector<double> source, cell &flow, 
        vector<int> wells)
{   
    boost::array<char, 3> border_flag = {0, 0, 0};
    double tmp = 0.;

    for(int index = 0; index < n; ++index)
    {
        tmp = 0.;

        if(!(is_well(index, wells)))
        {
            border_flag = check_border(index);

            //x-axis
            if(border_flag[0] != 'w')
                tmp += 
                    flow.x_left[index] * (conc_1[index] - conc_1[index - h_i]);

            if(border_flag[0] != 'e')
                tmp += 
                    flow.x_right[index] * (conc_1[index + h_i] - conc_1[index]);

            //y-axis
            if(border_flag[1] != 's')
                tmp += 
                    flow.y_left[index] * (conc_1[index] - conc_1[index - h_j]);

            if(border_flag[1] != 'n')
                tmp += 
                    flow.y_right[index] * (conc_1[index + h_j] - conc_1[index]);

            //z-axis
            if(border_flag[2] != 'b')
                tmp += 
                    flow.z_left[index] * (conc_1[index] - conc_1[index - h_k]);

            if(border_flag[2] != 0)
                tmp += 
                    flow.z_right[index] * (conc_1[index + h_k] - conc_1[index]);


            conc_2[index] = 
                (conc_1[index] + conc_2[index]) / 2. +
                ( h_t * 0.5 * time_un / porosity[index] ) * 
                ( -source[index] - tmp / (2. * length * h) );

            if(conc_2[index] < 0.)
                conc_2[index] = 0.;
        }
        else
            conc_1[index] = conc_2[index] = 0.;
    }

    return 0;
}

int flow_calc(vector<double> &pressure, cell &flow, 
        vector<double> &permeability, vector<int> wells)
{
    boost::array<char, 3> border_flag = {0, 0, 0};

    for(int index = 0; index < n; ++index)
    {
        border_flag = check_border(index);

        //x-axis flow
        if((border_flag[0] != 'w'))
            flow.x_left[index] = pressure[index] - pressure[index - h_i];
        else
            flow.x_left[index] = 0.;

        if((border_flag[0] != 'e'))
            flow.x_right[index] = pressure[index + h_i] - pressure[index];
        else
            flow.x_right[index] = 0.;

        //y-axis flow
        if((border_flag[1] != 's'))
            flow.y_left[index] = pressure[index] - pressure[index - h_j];
        else
            flow.y_left[index] = 0.;

        if((border_flag[1] != 'n'))
            flow.y_right[index] = pressure[index + h_j] - pressure[index];
        else
            flow.y_right[index] = 0.;

        //z-axis flow
        if((border_flag[2] != 'b'))
            flow.z_left[index] = pressure[index] - pressure[index - h_k];
        else
            flow.z_left[index] = 0.;

        if((border_flag[2] != 't'))
            flow.z_right[index] = pressure[index + h_k] - pressure[index];
        else
            flow.z_right[index] = 0.;

        flow.x_left[index] *= (-permeability[index] / (h * length));
        flow.y_left[index] *= (-permeability[index] / (h * length));
        flow.z_left[index] *= (-permeability[index] / (h * length));

        flow.x_right[index] *= (-permeability[index] / (h * length));
        flow.y_right[index] *= (-permeability[index] / (h * length));
        flow.z_right[index] *= (-permeability[index] / (h * length));
    }
    return 0;
}

int build_press_mat(vector<int> &col, vector<double> &val, vector<int> &ptr, 
        vector<double> &rhs, vector<double> &permeability, vector<int> &wells,
        vector<double> &source, vector<double> &dil_dt)
{    
    std::fill(rhs.begin(), rhs.end(), 0.);

    col.clear();  ptr.clear(); val.clear();  
    
    ptr.push_back(0);

    boost::array<char, 3> border_flag = {0, 0, 0};

    double k_f, k_b, j_f, j_b, i_f, i_b, cntr;

    for(int index = 0; index < n; ++index)
    {
        rhs[index] = h * h * (time_un * source[index] - dil_dt[index]);

        cntr = k_f = k_b = j_f = j_b = i_b = i_f = 0.;

        border_flag = check_border(index);

        if(border_flag[0] != 0)
        {
            /*if(border_flag[0] == 'w')
                i_f = 2. * get_pr_coef(index, index + h_i, permeability, wells);
            else
                i_b = 2. * get_pr_coef(index, index - h_i, permeability, wells);*/
            rhs[index] = p_top;
            cntr = 1.;
            k_f = k_b = j_f = j_b = i_b = i_f = 0.;
        }
        else
        {
            i_f = get_pr_coef(index, index + h_i, permeability, wells);
            i_b = get_pr_coef(index, index - h_i, permeability, wells);
        }               

        if(border_flag[1] != 0)
        {
            /*if(border_flag[1] == 's')
                j_f = 2. * get_pr_coef(index, index + h_j, permeability, wells);
            else
                j_b = 2. * get_pr_coef(index, index - h_j, permeability, wells);*/            
            rhs[index] = p_top;
            cntr = 1.;
            k_f = k_b = j_f = j_b = i_b = i_f = 0.;

        }
        else
        {
            j_f = get_pr_coef(index, index + h_j, permeability, wells);
            j_b = get_pr_coef(index, index - h_j, permeability, wells);
        }

        if(border_flag[2] != 0)
        {
            if(border_flag[2] == 'b')
                k_f = 2. * get_pr_coef(index, index + h_k, permeability, wells);
            else
                k_b = 2. * get_pr_coef(index, index - h_k, permeability, wells);
            /*rhs[index] = p_top;
            cntr = 1.;
            k_f = k_b = j_f = j_b = i_b = i_f = 0.;*/
        }
        else
        {
            k_f = get_pr_coef(index, index + h_k, permeability, wells);
            k_b = get_pr_coef(index, index - h_k, permeability, wells);
        }

        cntr += -(k_f + k_b + j_b + j_f + i_b + i_f);

        if(is_well(index, wells))
        {
            rhs[index] = p_bh;
            cntr = 1.;
            k_f = k_b = j_f = j_b = i_b = i_f = 0.;
        }

        if(is_well(index + h_j, wells))
        {
            rhs[index] -= j_f * p_bh;
            j_f = 0.;
        }

        if(is_well(index - h_j, wells))
        {
            rhs[index] -= j_b * p_bh;
            j_b = 0.;
        }

        if(is_well(index + h_i, wells))
        {
            rhs[index] -= i_f * p_bh;
            i_f = 0.;
        }

        if(is_well(index - h_i, wells))
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

    return 0;
}

double get_pr_coef(int index_1, int index_2, vector<double> &permeability, 
                        vector<int> &wells)
{
    int k_1 = index_1 / (n_x * n_y);
    int tmp_loc = index_1 - k_1 * n_x * n_y;
    int j_1 = tmp_loc / n_x;
    int i_1 = tmp_loc - j_1 * n_x;

    int k_2 = index_2 / (n_x * n_y);
    tmp_loc = index_2 - k_2 * n_x * n_y;
    int j_2 = tmp_loc / n_x;
    int i_2 = tmp_loc - j_2 * n_x;

    double perm = ( permeability[index_1] + permeability[index_2] ) / 2.;

    double tmp = perm * time_un / (length * length); 
 
    int dx = i_2 - i_1;
    int dy = j_2 - j_1;

    if(dx)
        tmp *= get_nu(double(i_1) + dx * 0.5, double(j_1), wells);
     
    if(dy)
        tmp *= get_nu(double(i_1), double(j_1) + dy * 0.5, wells);


    return -tmp;
}

boost::array<char, 3> check_border(int index)
{
    int k = index / (n_x*n_y);
    int tmp_loc = index - k*n_x*n_y;
    int j = tmp_loc/n_x;
    int i = tmp_loc - j*n_x;
    boost::array<char, 3> flag = {0, 0, 0};

    if(i == 0)
        flag[0] = 'w';
    
    if(i == n_x - 1)
        flag[0] = 'e';

    if(j == 0)
        flag[1] = 's';

    if(j == n_y - 1)
        flag[1] = 'n';

    if(k == 0)
        flag[2] = 'b';

    if(k == n_z - 1)
        flag[2] = 't';

    return flag;
}

bool is_border(int index)
{
    int k = index / (n_x*n_y);
    int tmp_loc = index - k*n_x*n_y;
    int j = tmp_loc/n_x;
    int i = tmp_loc - j*n_x;

    return ( (k == 0) || (k == n_z - 1) || (j == 0) || (j == n_y - 1) 
            || (i == 0) || (i == n_x - 1) );
}

bool is_well(int index, vector<int> &wells)
{
    return binary_search(wells.begin(), wells.end(), index);
}

int place_val(double value, int index, vector<int> &col, vector<double> &val)
{
    if(value)
    {
        val.push_back(value);
        col.push_back(index);
    }

    return 0;
}

int build_disp_mat(vector<int> &col,vector<double> &val, 
        vector<int> &ptr, vector<int> &wells)
{
    double c_1 = 1. + lame_1/lame_2;
    double c_2= 1. + c_1;          

    ptr.push_back(0);
    int point;

    for(int index = 0; index < n; ++index)
    {
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
    return 0;
}

int fill_disp_rhs(vector<double> &pressure, 
        vector<double> &rhs, vector<int> &wells)
{
    double tmp = 0.;

    for(int index = 0; index < n; ++index)
    {
        if( !(is_well(index, wells)) && !(is_border(index)) )
        {
            if(pressure[index] < pressure[index - h_i])
                tmp = (pressure[index] - pressure[index - h_i]);
            else
                tmp = (pressure[index + h_i] - pressure[index]);
            rhs[3*index + 0] = -h * tmp * length;

            if(pressure[index] < pressure[index - h_j])
                tmp = (pressure[index] - pressure[index - h_j]);
            else
                tmp = (pressure[index + h_j] - pressure[index]);
            rhs[3*index + 1] = -h * tmp * length;

            if(pressure[index] < pressure[index - h_k])
                tmp = (pressure[index] - pressure[index - h_k]);
            else
                tmp = (pressure[index + h_k] - pressure[index]);               
            rhs[3*index + 2] = -h * tmp * length;

        }
        else
            rhs[index] = 0.;
    }

    return 0;
}

int dil_calc(vector<double> &disp, vector<double> &dilatation, 
        vector<double> &dil_dt, vector<int> &wells)
{
    double tmp = 0.;

    for(int index = 0; index < n; ++index)
    {
        tmp = 0.;

        if( !(is_border(index)) && !(is_well(index, wells)) )
        {
            if(!(is_well(index + h_i, wells)) && !(is_well(index - h_i, wells)))
                tmp += disp[X(index + h_i)] - disp[X(index - h_i)];
            
            if(!(is_well(index + h_j, wells)) && !(is_well(index - h_j, wells)))
                tmp += disp[Y(index + h_j)] - disp[Y(index - h_j)];
    
            tmp += disp[Z(index + h_k)] - disp[Z(index - h_k)];

            tmp /= (2. * h);

            dil_dt[index] = ( tmp - dilatation[index] ) / h_t;

            dilatation[index] = tmp;
        }
        else
        {
            dilatation[index] = 0.;
            dil_dt[index] = 0.;
        }
    }

    return 0;
}


int por_calc(vector<double> &concentration, vector<double> &porosity, 
        vector<double> &source, cell &flow, 
        vector<double> &dil_dt, vector<int> &wells)
{
    for(int index = 0; index < n; ++index)
    {
        if( !(is_well(index, wells)) )
        {
            double vel_mod = 
                sqrt( 
                        pow(flow.x_left[index], 2.) + 
                        pow(flow.y_left[index], 2.) +
                        pow(flow.z_left[index], 2.) +
                        pow(flow.x_right[index], 2.) + 
                        pow(flow.y_right[index], 2.) +
                        pow(flow.z_right[index], 2.) 
                    );

            double tearoff = (vel_mod < q_0 ? 0. : alfa);

            source[index] = beta * concentration[index] 
                - tearoff * vel_mod / ( porosity[index] * ro_s ) 
                + tearoff * q_0 / ro_s;

            porosity[index] += 
                h_t * ( 
                        ( 1. - porosity[index] ) * dil_dt[index] 
                        - time_un * source[index] 
                      ); 

            if(porosity[index] > 0.7)
                porosity[index] = 0.7;
        }
        else
            porosity[index] = fi_0;
    }

    return 0;
}


int per_calc(vector<double> &porosity, vector<double> &permeability, 
        vector<int> &wells)
{
    for(int index = 0; index < n; ++index)
        if( !(is_well(index, wells)) )
        {
            permeability[index] *= pow(porosity[index], 3) 
                / pow(1. - porosity[index], 2);
        }

    return 0;
}


int add_well(int x, int y, vector<int> &wells)
{
    for(int k = 0; k < n_z; ++k)
        wells.push_back(x + y * n_x + k * n_x * n_y);
}


int drawLine(int x1, int y1, int x2, int y2, int z, std::vector<double> &matrix) 
{
    double frac = 100.0;
    //
    const int deltaX = abs(x2 - x1);
    const int deltaY = abs(y2 - y1);
    const int signX = x1 < x2 ? 1 : -1;
    const int signY = y1 < y2 ? 1 : -1;
    //
    int error = deltaX - deltaY;
    //
    if(matrix[x2 + z*n_x + y2*n_x*n_y]<(frac*k_0))
        matrix[x2 + z*n_x + y2*n_x*n_y] *= frac;

    while(x1 != x2 || y1 != y2) 
    {
        if(matrix[x1 + z*n_x + y1*n_x*n_y]<(frac*k_0))
            //if((x1<3)||(x1>5))
                matrix[x1 + z*n_x + y1*n_x*n_y] *= frac;

        const int error2 = error * 2;
        //
        if(error2 > -deltaY) {
            error -= deltaY;
            x1 += signX;
        }
        if(error2 < deltaX) {
            error += deltaX;
            y1 += signY;
        }
    }
    return 0;
}

inline int get_idx(int i, int j, int k)
{
    return i + j * n_x + k * n_x * n_y;
}

inline double edge_center(std::vector<double> &val, int index, int d1)
{
    return 0.5 * (val[index] + val[index + d1]);
}

inline double face_center(std::vector<double> &val, int index, int d1, int d2)
{

    double tmp  = val[index] 
                + val[index + d1]  
                + val[index + d2] 
                + val[index + d1 + d2];

    return 0.25 * tmp;
}

inline double cell_center(std::vector<double> &val, int index)
{
    double tmp  = val[index] 
                + val[index + h_i]  
                + val[index + h_j] 
                + val[index + h_k] 
                + val[index + h_i + h_j] 
                + val[index + h_i + h_k] 
                + val[index + h_j + h_k] 
                + val[index + h_i + h_j + h_k];

        return 0.125 * tmp;                                                               
}

std::vector<int> get_index(
        int i_start, 
        int i_end, 
        int j_start, 
        int j_end, 
        int k_start, 
        int k_end)
{
    std::vector<int> tmp;

    for(auto k = k_start; k < k_end; ++k)
        for(auto j = j_start; j < j_end; ++j)
            for(auto i = i_start; i < i_end; ++i)
                tmp.push_back(i + j * n_x + k * n_x * n_y);

    return tmp;
}

int get_flow(std::vector<double> &q, std::vector<double> &p,
        std::vector<double> &K)
{
    double qx, qy, qz;
    for(auto k = 1; k < n_z - 1; ++k)
        for(auto j = 1; j < n_y - 1; ++j)
            for(auto i = 1; i < n_x - 1; ++i)
            {
                auto idx = get_idx(i, j , k);   
                qx = undim * K[idx] * (p[idx + h_i] - p[idx - h_i]);
                qy = undim * K[idx] * (p[idx + h_j] - p[idx - h_j]);
                qx = undim * K[idx] * (p[idx + h_k] - p[idx - h_k]);
                q[idx] = sqrt(qx*qx + qy*qy + qz*qz);
            }

    return 0;    
}


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
    q1 = K_tmp * undim * (p2 - p1) / phi_tmp;

    p2 = edge_center(p, idx + d1, d2);
    p1 = edge_center(p, idx, d2);
    q2 = K_tmp * undim * (p2 - p1) / phi_tmp;

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

cell get_conc_cell(
        std::vector<double> &c, 
        int idx, 
        std::vector<double> &p,
        std::vector<double> &K,
        std::vector<double> &phi
        )
{
    cell tmp;

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
    tmp_x = undim * ( p2 - p1 ) * ( tmp.x_right[4] - tmp.x_left[4] );
 
    p2 = face_center(p, idx + h_j, h_i, h_k);
    p1 = face_center(p, idx, h_i, h_k);
    tmp_y = undim * ( p2 - p1 ) * ( tmp.y_right[4] - tmp.y_left[4] );

    p2 = face_center(p, idx + h_k, h_i, h_j);
    p1 = face_center(p, idx, h_i, h_j);
    tmp_z = undim * ( p2 - p1 ) * ( tmp.z_right[4] - tmp.z_left[4] );
   
    tmp.center = cell_center(c, idx) + K_tmp * h_t * (tmp_x + tmp_y + tmp_z) / (2. * h * phi_tmp);

    return tmp;
}

int get_c_vol(
        std::vector<cell> &c_vol,
        std::vector<double> &c, 
        std::vector<double> &p,
        std::vector<double> &K,
        std::vector<double> &phi
        )
{
    for(auto k = 1; k < n_z - 1; ++k)
        for(auto j = 1; j < n_y - 1; ++j)
            for(auto i = 1; i < n_x - 1; ++i)
            {
                auto idx = get_idx(i, j, k);

                c_vol[idx] = get_conc_cell(c, idx, p, K, phi);
            }

    return 0;
}

int get_c_flux(std::vector<cell> &c_vol, 
        std::vector<double> &flux, 
        std::vector<double> &p, 
        std::vector<double> &K, 
        std::vector<double> &phi, 
        int d1, 
        int d2,
        int d3
        )
{
    double tmp;
    double K_tmp, phi_tmp, p_tmp;

    for(auto k = 1; k < n_z - 1; ++k)
        for(auto j = 1; j < n_y - 1; ++j)
            for(auto i = 1; i < n_x - 1; ++i)
            {
                auto idx = get_idx(i, j, k);

                tmp = c_vol[idx].center 
                    + c_vol[idx - d1].center 
                    + c_vol[idx - d2].center
                    + c_vol[idx - d1 - d2].center;

                flux[idx] = 0.25 * tmp;
            }

    return 0;
}

int lax_wendroff_3d(
        std::vector<double> &c, 
        std::vector<cell> &c_vol, 
        std::vector<double> &p,
        std::vector<double> &K,
        std::vector<double> &phi,
        std::vector<int> &wells,
        std::vector<double> &q
        )
{
    double tmp;
    std::vector<double> F(n, 0.);
    std::vector<double> G(n, 0.);
    std::vector<double> V(n, 0.);
    double tmp_x, tmp_y, tmp_z;
    double p_tmp;
    double source;

    get_c_vol(c_vol, c, p, K, phi);
    
    get_c_flux(c_vol, F, p, K, phi, h_j, h_k, h_i);
    get_c_flux(c_vol, G, p, K, phi, h_i, h_k, h_j);
    get_c_flux(c_vol, V, p, K, phi, h_i, h_j, h_k);

    for(auto k = 1; k < n_z - 1; ++k)
        for(auto j = 1; j < n_y - 1; ++j)
            for(auto i = 1; i < n_x - 1; ++i)
            {
                auto idx = get_idx(i, j, k);

                p_tmp = ( p[idx + h_i] - p[idx - h_i] ) / 2.;
                tmp_x = ( F[idx] - F[idx - h_i] ) * undim * p_tmp * K[idx] / phi[idx];

                p_tmp = ( p[idx + h_j] - p[idx - h_j] ) / 2.;
                tmp_y = ( G[idx] - G[idx - h_j] ) * undim * p_tmp * K[idx] / phi[idx];

                p_tmp = ( p[idx + h_k] - p[idx - h_k] ) / 2.;
                tmp_z = ( V[idx] - V[idx - h_k] ) * undim * p_tmp * K[idx] / phi[idx];

                tmp = ( tmp_x + tmp_y + tmp_z ) * h_t / h;

                if(!is_well(idx, wells))
                {
                    c[idx] += tmp;

                    /*if(q[idx] / phi[idx] >= q_0)
                        source = beta * c[idx] - alfa * (q[idx] / phi[idx] - q_0);
                    else
                        source = beta * c[idx];
                    c[idx] -= h_t * source;*/

                    if((c[idx]) < 1e-9)
                        c[idx] = 0.;

                    if(c[idx] > 1.)
                        c[idx] = 1.;
                }
                else
                    c[idx] = 0.1;
            }

    return 0;
}
