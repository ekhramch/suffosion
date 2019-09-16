#include "utils.hpp"
#include "parameters.hpp"

Cell::Cell(int n = 5) : x_left(n, 0.), 
                        x_right(n, 0.),
                        y_left(n, 0.), 
                        y_right(n, 0.),
                        z_left(n, 0.), 
                        z_right(n, 0.),
                        center{0.}
{}

Unknowns::Unknowns(const int n) :   pressure(n, p_top),
                                    flow(n, 0),
                                    c_vol(n),
                                    source(n, 0),
                                    concentration(n, c_0),
                                    permeability(n, k_0),
                                    u_x(n, 0.), 
                                    u_y(n, 0.), 
                                    u_z(n, 0.),
                                    dilatation(n, 0.),
                                    dil_dt(n, 0.),
                                    porosity(n, fi_0),
                                    por_dt(n, 0.)
{
}

boost::array<char, 3> check_border(int index)
{
    int k = index / (n_x * n_y);
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

bool is_well(int i, int j, std::set<std::pair<int, int>> &wells)
{
    return wells.count(std::make_pair(i, j));
}

void place_val(double value, int index, std::vector<int> &col, std::vector<double> &val)
{
    if(value)
    {
        val.push_back(value);
        col.push_back(index);
    }
}