#ifndef PR_UTIL_H 
#define PR_UTIL_H

#include <boost/array.hpp>

const int n_x = 100; // amount of steps by x-axis
const int n_y = n_x; // amount of steps by x-axis
const int n_z = n_x; // amount of steps by x-axis
const int n_t = 1000; // amount of time steps
const int n = n_x * n_y * n_z; //size of vectors
const int h_i = 1; //for i+-1
const int h_j = n_x; //for j+-1
const int h_k = n_x * n_y; //for k+-1

const double h = 1. / double(n_x); // spatial step
const double h_t = 1. / double(n_t); // time step

const double time_un = 3600.; // time for undim, 1 hour, sec
const double length = 100.; // standart length of the layer, m
const double depth = 500.; // depth of layer, m

const double k_0 = 1e-12; // initial permeability, m^2=1darci
const double q_0 = 1.16*1e-5; // initial speed, 1 m per day in m/sec
const double fi_0 = 0.55; // initial porosity
const double c_0 = 0.; // initial concentration, kg/m^3

const double E = 7e6; // Young's modulus for soil, Pa
const double nu = 0.35; // Poisson's ratio, non-dim
const double lame_1 = E * nu / ( (1. - 2.*nu) * (1. + nu) ); // Lame first parameter(lambda), Pa
const double lame_2 = E / ( 2.*(1. + nu) ); // Lame second parameter(G), Pa

const double eta = 1e-3; // dynamic viscosity, Pa*s
const double ro_g = 1e3*9.81*length/lame_2; // density of water * grav acceleration(Earth), undimensioned
const double ro_s = 1.; // density of the solid phase, kg/m^3
const double alfa = 1.; // parameter of the process, kg/m^4
const double beta = 0.1; // parameter of the process, 1/sec
const double gamma_1 = 0.1; // parameter, m^2/s
const double d = 1e-4; // diameter of the particles, m
const double T = 10.; //tortuosity

const double p_top = 0.1*1e6 + 0.1*1e6*depth/10.; // initial pressure upper: atmospheric + 0.1 MPa for 10 m 
const double p_bot = p_top; //minus hydrostatic
//const double p_bot = 0.1*1e6 + 0.1*1e6*(depth+length)/10.; // initial pressure bottom: atmospheric + 0.1 MPa for 10 m 
const double p_bh = 2. * p_top; // pressure on the borehole
const double r_c = 0.1 / length; // radius of borehole, undimensioned

struct cell
{
    std::vector<double> x_left;
    std::vector<double> x_right;
    
    std::vector<double> y_left;
    std::vector<double> y_right;
    
    std::vector<double> z_left;   
    std::vector<double> z_right;

    double center;

    cell(int n = 5) : x_left(n, 0.), 
                      x_right(n, 0.),
                      y_left(n, 0.), 
                      y_right(n, 0.),
                      z_left(n, 0.), 
                      z_right(n, 0.),
                      center{0.}
    {}
};

int wrt_vtk(std::vector<double> &arr, const std::string &filename);

double get_nu(double i, double j, std::vector<int> wells);

int flow_calc(std::vector<double> &pressure, cell &flow,
        std::vector<double> &permeability, std::vector<int> wells);

int conc_calc(std::vector<double> &conc_1, std::vector<double> &conc_2,
        std::vector<double> porosity, std::vector<double> source, 
        cell &flow, std::vector<int> wells);

int build_press_mat(std::vector<int> &col, std::vector<double> &val, 
        std::vector<int> &ptr, std::vector<double> &rhs, 
        std::vector<double> &permeability, std::vector<int> &wells,
        std::vector<double> &source, std::vector<double> &dil_dt);

int drawLine(int x1, int y1, int x2, int y2, int z,std::vector<double> &matrix); 

double get_perm(int index_1, int index_2, std::vector<double> &permeability);

double get_pr_coef(int index_1, int index_2, std::vector<double> &permeability,
        std::vector<int> &wells);

boost::array<char, 3> check_border(int index);

bool is_border(int index);

bool is_well(int index, std::vector<int> &wells);

int place_val(double value, int index, 
        std::vector<int> &col, std::vector<double> &val);

int build_disp_mat(std::vector<int> &col, std::vector<double> &val, 
        std::vector<int> &ptr, std::vector<int> &wells);

int fill_disp_rhs(std::vector<double> &pressure, std::vector<double> &rhs,
        std::vector<int> &wells);

int dil_calc(std::vector<double> &disp, std::vector<double> &dilatation, 
        std::vector<double> &dil_dt, std::vector<int> &wells);

int por_calc(std::vector<double> &concentration, 
        std::vector<double> &porosity, std::vector<double> &source,
        cell &flow, std::vector<double> &dil_dt, std::vector<int> &wells);

int per_calc(std::vector<double> &porosity, std::vector<double> &permeability,
        std::vector<int> &wells);

int add_well(int x, int y, std::vector<int> &wells);

int get_flow(std::vector<cell> &q, std::vector<double> &p, 
        std::vector<double> &perm);

int lax_wendroff_3d(std::vector<double> &c, std::vector<cell> &c_vol, 
        std::vector<cell> &q);

#endif
