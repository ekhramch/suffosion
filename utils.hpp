#pragma once
#include <vector>
#include <utility>
#include <set>
#include "boost/multi_array.hpp"

//structure for flux calc
struct Cell
{
    std::vector<double> x_left;
    std::vector<double> x_right;
    
    std::vector<double> y_left;
    std::vector<double> y_right;
    
    std::vector<double> z_left;   
    std::vector<double> z_right;

    double center;

    Cell(int n);
};

//wrapper for vectors of unknowns
struct Unknowns
{
    std::vector<double> pressure; //pressure in water
    std::vector<double> flow; //abs of velocities = |q|
    std::vector<Cell> c_vol; //volume concentration
    std::vector<double> source;//source of solid phase
    std::vector<double> concentration; //concentration of particles
    std::vector<double> permeability; //distribution of permeability
    std::vector<double> u_x, u_y, u_z; //displacements
    std::vector<double> dilatation; //trace of strain tensor
    std::vector<double> dil_dt; //dilatation time derivative
    std::vector<double> porosity;
    std::vector<double> por_dt; //porosity time derivative
    std::set<std::pair<int, int>> wells; //coordinates of wells in 2d

    Unknowns(const int n);
};

boost::array<char, 3> check_border(int index);

bool is_border(int index);

bool is_well(int i, int j, std::set<std::pair<int, int>> &wells);

void place_val(double value, int index, std::vector<int> &col, std::vector<double> &val);