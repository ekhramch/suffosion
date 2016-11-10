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

    f << "POINT_DATA " << (n_x) * (n_y) * n_z << endl;

    f << "SCALARS data double" << endl;

    f << "LOOKUP_TABLE default" << endl;

    for(int k=0; k<n_z; k++)
    {
	for(int j=0; j<n_y; j++)
	    for(int i=0; i<n_x; i++)
		f <<  arr[i + j*n_x + k*n_x*n_y] << " ";

	f << endl;
    }

    return 0;
}

double get_nu(double i, double j, vector<point> boreholes)
{
    auto sum = 0.;

    for(auto& bh : boreholes)
    {
        auto i_b = bh.x;
        auto j_b = bh.y;
        auto r_2 = (double(i_b) - i)*(double(i_b) - i) + (double(j_b) - j)*(double(j_b) - j);

        if(fabs(r_2 - 0.25)<0.0001)
            sum += 3.14/(2.0*log((h/r_c)));
        else
            if(fabs(r_2 - 1.25)<0.0001)
                sum += 2.0*atan(0.5)/log(2.0);
            else
                if(fabs(r_2 - 2.25)<0.0001)
                    sum += 2.0*atan(1.0/3.0)/log(2.0);
    }

    return (sum ? sum : 1.);
}

int get_q(vector<double> &pressure, vec_3d &velocity, vector<double> &permeability, vector<point> boreholes)
{
    for(auto index = 0; index < n; ++index)
    {
        auto k = index / (n_x * n_y);
        auto tmp_loc = index - (k * n_x * n_y);
        auto j = tmp_loc / n_x;
        auto i = tmp_loc - (j * n_x);

        if ((i == 0) || (i == n_x))
            velocity.x[index] = 0.;
        else
        {
            if (pressure[index] > pressure[index - 1])
                velocity.x[index] = (pressure[index + 1] - pressure[index]) / h;
            else
                velocity.x[index] = (pressure[index] - pressure[index - 1]) / h;
        }

        if ((j == 0) || (j == n_y))
            velocity.y[index] = 0.;
        else
        {
            if (pressure[index] > pressure[index - n_x])
                velocity.y[index] = (pressure[index + n_x] - pressure[index]) / h;
            else
                velocity.y[index] = (pressure[index] - pressure[index - n_x]) / h;
        }

        if ((k == 0) || (k == n_z))
            velocity.z[index] = 0.;
        else
        {
            if (pressure[index] > pressure[index - n_x*n_y])
                velocity.z[index] = (pressure[index + n_x*n_y] - pressure[index]) / h;
            else
                velocity.z[index] = (pressure[index] - pressure[index - n_x*n_y]) / h;
        }

        velocity.x[index] *= -permeability[index] / length;
        velocity.y[index] *= -permeability[index] / length;
        velocity.z[index] *= -permeability[index] / length;
    }

    return 0;
}

int set_bc(vector<double> &q)
{

    for(int k=0, n=1; k<n_z; k+=n_z-1, n-=2)
#pragma omp parallel for
        for(int j = 0; j < n_y; ++j)
            for (int i = 0; i < n_x; ++i)
                q[i + j*n_x + k*n_x*n_y] = q[i + j*n_x + (k+n)*n_x*n_y];

    for(int j=0, n=1; j<n_y; j+=n_y-1, n-=2)
#pragma omp parallel for
        for(int k = 1; k < n_z-1; ++k)
            for (int i = 0; i < n_x; ++i)
                q[i + j*n_x + k*n_x*n_y] = q[i + (j+n)*n_x + k*n_x*n_y];


    for(int i=0, n=1; i<n_x; i+=n_x-1, n-=2)
#pragma omp parallel for
        for(int k = 1; k < n_z-1; ++k)
            for(int j = 1; j < n_y-1; ++j)
                q[i + j*n_x + k*n_x*n_y] = q[(i+n) + j*n_x + k*n_x*n_y];
    return 0;
}


int build_press_mat(vector<int> &col, vector<double> &val, vector<int> &ptr, 
        vector<double> &rhs, vector<double> &permeability, vector<point> &boreholes)
{

    ptr.push_back(0);

    //Build matrix
    for(auto index=0; index < n; ++index)
    {
        auto k = index/(n_x*n_y);
        auto tmp_loc = index - k*n_x*n_y;
        auto j = tmp_loc/n_x;
        auto i = tmp_loc - j*n_x;

        point tmp_bh(i,j);
        bool is_bh = ( (!boreholes.empty()) && (find(boreholes.begin(), boreholes.end(), tmp_bh) != boreholes.end()) );

        auto undim = hour / (length * length); //undimensional coeff for d_2_p/d_xi_2, pressure still in pascals


        auto k_f = -undim * ( permeability[index] + permeability[index - n_x*n_y] )/2.;
        auto k_b = -undim * ( permeability[index] + permeability[index + n_x*n_y] )/2.;
        auto j_f = -undim * get_nu(double(i), double(j + 0.5), boreholes)*( permeability[index] + permeability[index + n_x] )/2.;
        auto j_b = -undim * get_nu(double(i), double(j - 0.5), boreholes)*( permeability[index] + permeability[index - n_x] )/2.;
        auto i_f = -undim * get_nu(double(i + 0.5), double(j), boreholes)*( permeability[index] + permeability[index + 1] )/2.;
        auto i_b = -undim * get_nu(double(i - 0.5), double(j), boreholes)*( permeability[index] + permeability[index - 1] )/2.;
        auto cntr = -(k_f + k_b + j_f + j_b + i_f + i_b);

        if(k==0)
        {
            k_f = k_b = j_f = j_b = i_f = i_b = 0;
            cntr = 1.;
            rhs[index] = p_up;
        }

        if(k==n_z-1)
        {
            k_f = k_b = j_f = j_b = i_f = i_b = 0;
            cntr = 1.;
            rhs[index] = p_bot;
        }

        if( (j==0) && j_b )
        {
            j_f *= 2;
            j_b = 0;
        }

        if( (j==n_y-1) && j_f )
        {
            j_f = 0;
            j_b *= 2;
        }

        if( (i==0) && i_b )
        {
            i_f *= 2;
            i_b = 0;
        }

        if( (i==n_x-1) && i_f )
        {
            i_f = 0;
            i_b *= 2;
        }

        if(is_bh)
        {
            k_f = k_b = j_f = j_b = i_f = i_b = 0;
            cntr = 1.;
            rhs[index] = p_bh;
        }

        if(k_b)
        {
            val.push_back(k_b);
            col.push_back(index - n_x*n_y);
        }

        tmp_bh.y = j - 1;
        is_bh = ( (!boreholes.empty()) && (find(boreholes.begin(), boreholes.end(), tmp_bh) != boreholes.end()) );

        if(j_b)
        {
            if(is_bh)
                rhs[index] -= j_b*p_bh;
            else
            {
                val.push_back(j_b);
                col.push_back(index - n_x);
            }
        }

        tmp_bh.y = j;
        tmp_bh.x = i - 1;
        is_bh = ( (!boreholes.empty()) && (find(boreholes.begin(), boreholes.end(), tmp_bh) != boreholes.end()) );

        if(i_b)
        {
            if(is_bh)
                rhs[index] -= i_b*p_bh;
            else
            {
                val.push_back(i_b);
                col.push_back(index - 1);
            }
        }

        val.push_back(cntr);
        col.push_back(index);

        tmp_bh.x = i + 1;
        is_bh = ( (!boreholes.empty()) && (find(boreholes.begin(), boreholes.end(), tmp_bh) != boreholes.end()) );

        if(i_f)
        {        
            if(is_bh)
                rhs[index] -= i_f*p_bh;
            else
            {
                val.push_back(i_f);
                col.push_back(index + 1);
            }
        }

        tmp_bh.y = j + 1;
        tmp_bh.x = i;
        is_bh = ( (!boreholes.empty()) && (find(boreholes.begin(), boreholes.end(), tmp_bh) != boreholes.end()) );

        if(j_f)
        {            
            if(is_bh)
                rhs[index] -= j_f*p_bh;
            else
            {
                val.push_back(j_f);
                col.push_back(index + n_x);
            }
        }

        if(k_f)
        {
            val.push_back(k_f);
            col.push_back(index + n_x*n_y);
        }

        ptr.push_back(col.size());
    }

    return 0;
}

int build_u_mat(vector<int> &col, vector<double> &val, vector<int> &ptr, char dir)
{
    double d_x = -1., d_y = -1., d_z = -1.;

    switch(dir)
    {    
        case 'x':
            d_x = -2 - lame_1/lame_2;
            break;
        case 'y':
            d_y = -2 - lame_1/lame_2;
            break;
        case 'z':
            d_z = -2 - lame_1/lame_2;
            break;
        default:
            return -1;
    }

    ptr.push_back(0);

    //Build matrix
    for(auto index = 0; index < n; ++index)
    {
        auto k = index / (n_x*n_y);
        auto tmp_loc = index - k*n_x*n_y;
        auto j = tmp_loc/n_x;
        auto i = tmp_loc - j*n_x;


        if ( ( k > 0 && k < ( n_z - 1 ) ) && ( j > 0 && ( j < n_y - 1 ) ) && ( i > 0 && i< ( n_x - 1 ) ) )
        {
            val.push_back(d_z);
            col.push_back(index - n_x*n_y);

            val.push_back(d_y);
            col.push_back(index - n_x);

            val.push_back(d_x);
            col.push_back(index - 1);

            val.push_back(8. + 2.*lame_1/lame_2);
            col.push_back(index);

            val.push_back(d_x);
            col.push_back(index + 1);

            val.push_back(d_y);
            col.push_back(index + n_x);

            val.push_back(d_z);
            col.push_back(index + n_x*n_y);
        }
        else
        {
            val.push_back(1.);
            col.push_back(index);
        }

        ptr.push_back(col.size());
    }

    return 0;
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

//mixed second order derivative; direction is determined by dir
double sec_ord_mx(vector<double> &u, int i, int j, int k, const int dir)
{
    double tmp=0.;
    auto index = i + j*n_x + k*n_x*n_y;

    //-45 for dxdy; 45 for dxdz; 135 for dydz

    switch(dir)
    {    
        case -45:
            //tmp = u[index + n_x] - u[index - 1 + n_x] - u[index] + u[index - 1];
            tmp = u[index + 1 + n_x] - u[index + 1 - n_x] - u[index - 1 + n_x] + u[index - 1 - n_x];
            break;
        case 45:
            //tmp = u[index + n_x*n_y] - u[index - 1 + n_x*n_y] - u[index] + u[index - 1];
            tmp = u[index + 1 + n_x*n_y] - u[index + 1 - n_x*n_y] - u[index - 1 + n_x*n_y] + u[index - 1 - n_x*n_y];
            break;
        case 135:
            //tmp = u[index + n_x*n_y] - u[index - n_x + n_x*n_y] - u[index] + u[index - n_x];
            tmp = u[index + n_x + n_x*n_y] - u[index + n_x - n_x*n_y] - u[index - n_x + n_x*n_y] + u[index - n_x - n_x*n_y];
            break;
        default:
            return -1;
    }

    return tmp/(4. * h * h);
}

