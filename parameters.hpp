#pragma once

//geometry
extern const double depth; //depth of the layer, m
extern const int L_x; //x-dimension of the layer, m
extern const int L_y; //y-dimension of the layer, m
extern const int L_z; //z-dimension of the layer, m
extern const int L_t; //typical time of 1 hour, s
extern const int n_x; // amount of steps by x-axis
extern const int n_y; // amount of steps by x-axis
extern const int n_z; // amount of steps by x-axis
extern const int n_t; // amount of time steps
extern const int n;// * n_y * n_z; //size of vectors
extern const int h_i; //for i+-1
extern const int h_j; //for j+-1
extern const int h_k; //for k+-1
extern const double h; // spatial step
extern const double h_t; // time step
extern const double depth; // depth of layer, m


//physical properties
extern const double E; // Young's modulus for soil, Pa
extern const double nu; // Poisson's ratio, non-dim
extern const double lame_1; // Lame first parameter(lambda), Pa
extern const double lame_2; // Lame second parameter(G), Pa
extern const double undim; //undimensional constant for Lame system
extern const double eta; // dynamic viscosity, Pa*s
extern const double ro_g; // density of water * grav acceleration(Earth), undimensioned
extern const double ro_s; // density of the solid phase, kg/m^3
extern const double alfa; // parameter of the process, kg/m^4
extern const double beta; // parameter of the process, 1/sec
extern const double b; // parameter of the process
extern const double gamma_1; // parameter, m^2/s
extern const double d; // diameter of the particles, m
extern const double T; //tortuosity
extern const double c_0; // initial concentration, kg/m^3
extern const double fi_0; // initial porosity
extern const double fi_max; // maximum porosity
extern const double grain_size; //between very fine sand and slit * [micrometer to mm]
extern const double unity_factor; //alpha (mdarcy/mm^2) from wiki * [mdarcy to darcy]
extern const double darcy_coef; // darcy to m^2
extern const double koz_car_coef; //constants of Cozeny-Karman law brought together
extern const double poros_coef; // porosity function
extern const double k_0; // initial permeability
extern const double q_0; // velocity of tearoff
extern const double p_top; // initial pressure upper: atmospheric + 0.1 MPa for 10 m 
extern const double p_bot; //minus hydrostatic
extern const double p_bh; // pressure on the borehole
extern const double r_c; // radius of borehole, undimensioned