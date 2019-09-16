#include "parameters.hpp"

//geometry
const double depth = 100.; //depth of the layer, m
const int L_x = 100; //x-dimension of the layer, m
const int L_y = 100; //y-dimension of the layer, m
const int L_z = 50; //z-dimension of the layer, m
const int L_t = 3600; //typical time of 1 hour, s
const int n_x = 50; // amount of steps by x-axis
const int n_y = 100; // amount of steps by x-axis
const int n_z = 100; // amount of steps by x-axis
const int n_t = 1000; // amount of time steps
const int n = n_x * n_y * n_z; //size of vectors
const int h_i = 1; //for i+-1
const int h_j = n_x; //for j+-1
const int h_k = n_x * n_y; //for k+-1
const double h = double(L_x) / double(n_x); // spatial step
const double h_t = double(L_t) / double(n_t); // time step
const double depth = 500.; // depth of layer, m

//physical properties
const double E = 7e6; // Young's modulus for soil, Pa
const double nu = 0.35; // Poisson's ratio, non-dim
const double lame_1 = E * nu / ( (1. - 2.*nu) * (1. + nu) ); // Lame first parameter(lambda), Pa
const double lame_2 = E / ( 2.*(1. + nu) ); // Lame second parameter(G), Pa
const double undim = lame_2 * time_un / ( length * length ); //undimensional constant for Lame system
const double eta = 1e-3; // dynamic viscosity, Pa*s
const double ro_g = 1e3*9.81*length; // density of water * grav acceleration(Earth), undimensioned
const double ro_s = 1.; // density of the solid phase, kg/m^3
const double alfa = 100.; // parameter of the process, kg/m^4
const double beta = 1e7; // parameter of the process, 1/sec
const double b = 0.01; // parameter of the process
const double gamma_1 = 0.1; // parameter, m^2/s
const double d = 1e-4; // diameter of the particles, m
const double T = 10.; //tortuosity
const double c_0 = 0.; // initial concentration, kg/m^3
const double fi_0 = 0.5; // initial porosity
const double fi_max = 0.75; // maximum porosity
const double grain_size = 50. * 1e-3; //between very fine sand and slit * [micrometer to mm]
const double unity_factor = ( (0.8 * 1e6) / 1.0135 ) * 1e-3; //alpha (mdarcy/mm^2) from wiki * [mdarcy to darcy]
const double darcy_coef = 9.869233e-13; // darcy to m^2
const double koz_car_coef = darcy_coef * grain_size * grain_size * unity_factor / eta ; //constants of Cozeny-Karman law brought together
const double poros_coef = fi_0 * fi_0 * fi_0 / ( ( 1. - fi_0 ) * ( 1. - fi_0 ) ); // porosity function
const double k_0 = koz_car_coef * poros_coef; // initial permeability
const double q_0 = 7e-5; // velocity of tearoff
const double p_top = (0.1*1e6 + 0.1*1e6*depth/10.) / lame_2; // initial pressure upper: atmospheric + 0.1 MPa for 10 m 
const double p_bot = p_top; //minus hydrostatic
const double p_bh = 2. * p_top; // pressure on the borehole
const double r_c = 0.1 / length; // radius of borehole, undimensioned