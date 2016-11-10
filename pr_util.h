#ifndef PR_UTIL_H 
#define PR_UTIL_H

const unsigned long n_x = 50; // amount of steps by x-axis
const unsigned long n_y = n_x; // amount of steps by x-axis
const unsigned long n_z = n_y; // amount of steps by x-axis
const unsigned long n_t = 1000; // amount of time steps
const unsigned long n = n_x*n_y*n_z; //size of vectors

const double h = 1./n_x; // spatial step
const double h_t = 1/n_t; // time step

const double hour = 3600.; // time for undim, 1 hour, sec
const double length = 100.; // standart length of the layer, m
const double depth = 500; // depth of layer, m

const double k_0 = 1e-12; // initial permeability, m^2=1darci
const double q_0 = 1.16*1e-5; // initial speed, 1 m per day in m/sec
const double fi_0 = 0.55; // initial porosity
const double c_0 = 1; // initial concentration, kg/m^3

const double E = 7e6; // Young's modulus for soil, Pa
const double nu = 0.35; // Poisson's ratio, non-dim
const double lame_1 = E * nu / ( (1. - 2.*nu) * (1. + nu) ); // Lame first parameter(lambda), Pa
const double lame_2 = E / ( 2.*(1. + nu) ); // Lame second parameter(G), Pa

const double eta = 1e-3; // dynamic viscosity, Pa*s
const double ro_g = 1e3*9.81*length/lame_2; // density of water * grav acceleration(Earth), undimensioned
const double ro_s = 1.; // density of the solid phase, kg/m^3
const double alfa = 0.1; // parameter of the process, kg/m^4
const double beta = 0.1; // parameter of the process, 1/sec
const double gamma_1 = 0.1; // parameter, m^2/s
const double d = 1e-4; // diameter of the particles, m
const double T = 10.; //tortuosity

const double p_up = 0.1*1e6 + 0.1*1e6*depth/10.; // initial pressure upper: atmospheric + 0.1 MPa for 10 m, undimensioned 
const double p_bot = (0.1*1e6 + 0.1*1e6*(depth+length)/10.); // initial pressure bottom: atmospheric + 0.1 MPa for 10 m, undimensioned 
const double p_bh = 2*p_bot; // pressure on the borehole
const double r_c = 0.1/length; // radius of borehole, undimensioned

struct point
{
    int x;
    int y;

    point(int a=0, int b=0) : x(a), y(b) {}

    public:
    inline bool operator==(const point &rhs)
    {
        return ( (x == rhs.x) && (y == rhs.y) );
    }

};

struct vec_3d
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    vec_3d(int n=1) : x(n,0.), y(n,0.), z(n,0.) {}

};

int wrt_vtk(std::vector<double> &arr, const std::string filename);

double get_nu(double i, double j, std::vector<point> boreholes);

int get_q(std::vector<double> &pressure, vec_3d &velocity, std::vector<double> &permeability, std::vector<point> boreholes);

int set_bc(std::vector<double> &q);

bool is_borehole(std::vector<point> &boreholes, const int i, const int j);

int build_press_mat(std::vector<int> &col, std::vector<double> &val, 
        std::vector<int> &ptr, std::vector<double> &rhs, 
        std::vector<double> &permeability, std::vector<point> &boreholes);

int build_u_mat(std::vector<int> &col, std::vector<double> &val, 
        std::vector<int> &ptr, char dir);

int drawLine(int x1, int y1, int x2, int y2, int z, std::vector<double> &matrix); 

double sec_ord_mx(std::vector<double> &arr, int i, int j, int k, const int dir);

#endif
