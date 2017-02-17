#include <vector>
#include <map>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <string>
#include "pr_util.h"
#include <sstream>
#include <assert.h>
#include "omp.h"

#define X(i)(3*(i) + 0)
#define Y(i)(3*(i) + 1)
#define Z(i)(3*(i) + 2)

using namespace std;

int wrt_vtk(vector<double> &arr, const string filename)
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

    for(auto index = 0; index < wells.size(); index += n_z)
    {
        auto z = index / (n_x * n_y);
        auto tmp_loc = index - z * n_x * n_y;
        auto y = tmp_loc / n_x;
        auto x = tmp_loc - y * n_x;

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
        vector<double> porosity, vector<double> source, 
        vec_3d &velocity, vector<int> wells)
{   
    array<char, 3> border_flag = {0, 0, 0};
    double tmp[3] = {0., 0., 0.};

    for(auto index = 0; index < n; ++index)
    {
        tmp[0] = tmp[1] = tmp[2] = 0.;

        int s_x = ( velocity.x[index] > 0 ? 0 : 1);
        int s_y = ( velocity.y[index] > 0 ? 0 : 1);
        int s_z = ( velocity.z[index] > 0 ? 0 : 1);

        border_flag = check_border(index);

        if(border_flag[0] == 0)
            tmp[0] = velocity.x[index] * 
                (conc_1[index + s_x * h_i] - conc_1[index - (!s_x) * h_i]);

        if(border_flag[1] == 0)
            tmp[1] = velocity.y[index] * 
                (conc_1[index + s_y * h_j] - conc_1[index - (!s_y) * h_j]);


        if(border_flag[2] == 0)
            tmp[2] = velocity.z[index] * 
                (conc_1[index + s_z * h_k] - conc_1[index - (!s_z) * h_k]);


        conc_2[index] = 
            (conc_1[index] + conc_2[index]) / 2. +
            ( h_t * 0.5 * time_un / porosity[index] ) * 
            ( -source[index] - (tmp[0] + tmp[1] + tmp[2]) / (length * h) );
    }

    return 0;
}

int vel_calc(vector<double> &pressure, 
        vec_3d &velocity, vector<double> &permeability, vector<int> wells)
{
    array<char, 3> border_flag = {0, 0, 0};

    for(auto index = 0; index < n; ++index)
    {
        border_flag = check_border(index);

        if(border_flag[0] == 0)
            if(pressure[index] > pressure[index - h_i])
                velocity.x[index] = (pressure[index + h_i] - pressure[index]);
            else
                velocity.x[index] = (pressure[index] - pressure[index - h_i]);

        if(border_flag[1] == 0)
            if(pressure[index] > pressure[index - h_j])
                velocity.y[index] = (pressure[index + h_j] - pressure[index]);
            else
                velocity.y[index] = (pressure[index] - pressure[index - h_j]);

        if(border_flag[2] == 0)
            if(pressure[index] > pressure[index - h_k])
                velocity.z[index] = (pressure[index + h_k] - pressure[index]);
            else
                velocity.z[index] = (pressure[index] - pressure[index - h_k]);


        velocity.x[index] *= -permeability[index] / ( h * length );

        velocity.y[index] *= -permeability[index] / ( h * length );

        velocity.z[index] *= -permeability[index] / ( h * length );
    }

    return 0;
}

int build_press_mat(vector<int> &col, vector<double> &val, vector<int> &ptr, 
        vector<double> &rhs, vector<double> &permeability, vector<int> &wells)
{    
    std::fill(rhs.begin(), rhs.end(), 0.);
    col.clear();  ptr.clear(); val.clear();  
    
    ptr.push_back(0);

    array<char, 3> border_flag = {0, 0, 0};

    double k_f, k_b, j_f, j_b, i_f, i_b, cntr;

    for(auto index = 0; index < n; ++index)
    {    
        cntr = k_f = k_b = j_f = j_b = i_b = i_f = 0.;

        border_flag = check_border(index);

        if(border_flag[0] != 0)
        {
            if(border_flag[0] == 'w')
                i_f = 2. * get_pr_coef(index, index + h_i, permeability, wells);
            else
                i_b = 2. * get_pr_coef(index, index - h_i, permeability, wells);
        }
        else
        {
            i_f = get_pr_coef(index, index + h_i, permeability, wells);
            i_b = get_pr_coef(index, index - h_i, permeability, wells);
        }               

        if(border_flag[1] != 0)
        {
            if(border_flag[1] == 's')
                j_f = 2. * get_pr_coef(index, index + h_j, permeability, wells);
            else
                j_b = 2. * get_pr_coef(index, index - h_j, permeability, wells);
        }
        else
        {
            j_f = get_pr_coef(index, index + h_j, permeability, wells);
            j_b = get_pr_coef(index, index - h_j, permeability, wells);
        }

        if(border_flag[2] != 0)
        {
            rhs[index] = p_top;
            cntr = 1.;
            k_f = k_b = j_f = j_b = i_b = i_f = 0.;
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

    if(dx == 1 || dx == -1)
        tmp *= get_nu(double(i_1) + dx * 0.5, double(j_1), wells);
     else if(dy == 1 || dy == -1)
        tmp *= get_nu(double(i_1), double(j_1) + dy * 0.5, wells);

    return -tmp;
}

array<char, 3> check_border(int index)
{
    auto k = index / (n_x*n_y);
    auto tmp_loc = index - k*n_x*n_y;
    auto j = tmp_loc/n_x;
    auto i = tmp_loc - j*n_x;

    array<char, 3> flag = {0, 0, 0};

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
    auto k = index / (n_x*n_y);
    auto tmp_loc = index - k*n_x*n_y;
    auto j = tmp_loc/n_x;
    auto i = tmp_loc - j*n_x;

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

    for(auto index = 0; index < n; ++index)
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
            if(!is_well(point, wells))
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
            if(!is_well(point, wells))
            {
                //u_y
                place_val(-c_1, Y(point), col, val);         
                //u_z
                place_val(-c_2, Z(point), col, val);         
            }

            ptr.push_back(col.size());
        }
    }
    return 0;
}

int fill_disp_rhs(vector<double> &pressure, 
        vector<double> &rhs, vector<int> &wells)
{
    double tmp = 0.;

    for(auto index = 0; index < n; ++index)
    {
        if( !(is_well(index, wells)) && !(is_border(index)) )
        {
            if(pressure[index] < pressure[index - h_i])
                tmp = (pressure[index] - pressure[index - h_i]);
            else
                tmp = (pressure[index + h_i] - pressure[index]);
            rhs[3*index + 0] = -h * tmp * length / lame_2;

            if(pressure[index] < pressure[index - h_j])
                tmp = (pressure[index] - pressure[index - h_j]);
            else
                tmp = (pressure[index + h_j] - pressure[index]);
            rhs[3*index + 1] = -h * tmp * length / lame_2;

            if(pressure[index] < pressure[index - h_k])
                tmp = (pressure[index] - pressure[index - h_k]);
            else
                tmp = (pressure[index + h_k] - pressure[index]);               
            rhs[3*index + 2] = -h * tmp * length / lame_2;

        }
    }

    return 0;
}

int dil_calc(vector<double> &disp, 
        vector<double> &dilatation, vector<double> &dil_dt)
{
    double tmp = 0.;

    for(auto index = 0; index < n; ++index)
    {
        if( !(is_border(index)) )
        {
            tmp  = disp[X(index + h_i)] - disp[X(index - h_i)];
            tmp += disp[Y(index + h_j)] - disp[Y(index - h_j)];
            tmp += disp[Z(index + h_k)] - disp[Z(index - h_k)];
            tmp /= (2. * h);

            dil_dt[index] = ( tmp - dilatation[index] ) / h_t;

            dilatation[index] = tmp;
        }
    }

    return 0;
}

int por_calc(vector<double> &concentration, vector<double> porosity, 
        vector<double> source, vec_3d &velocity, vector<double> dil_dt)
{
    for(auto index = 0; index < n; ++index)
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
    }

    return 0;
}

int per_calc(vector<double> &porosity, vector<double> permeability)
{
    for(auto index = 0; index < n; ++index)
    {
        double s = 6 * ( 1 - porosity[index] ) / d;

        permeability[index] = pow(porosity[index], 3)/(T * T * s * s * eta);
    }

    return 0;
}

int add_well(int x, int y, vector<int> &wells)
{
    for(auto k = 0; k < n_z; ++k)
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
