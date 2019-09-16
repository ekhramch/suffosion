#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cassert>

#include "saver.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "pressure.hpp"
#include "stress.hpp"
#include "chemistry.hpp"

/*AMGCL-VexCL include section*/
#include <vexcl/vexcl.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/io/mm.hpp>
/*end of section*/

namespace amgcl
{
    profiler<> prof("project suffosion");
}
using amgcl::prof;

int main(int argc, char *argv[])
{
    // Initialize VexCL context.
    vex::Context ctx( vex::Filter::Env && vex::Filter::DoublePrecision );
    if (!ctx)
    {
        std::cerr << "no GPUs" << std::endl;
        return 1;
    }
    std::cout << ctx << std::endl;

	//unknowns of the process
    Unknowns state(n);

	//hdf5 writer
    std::map<std::string, double*> save_data = 
    {
        {"pressure", state.pressure.data()},
        {"dilatation", state.dilatation.data()},
        {"porosity", state.porosity.data()},
        {"flow", state.flow.data()},
        {"permeability", state.permeability.data()},
        {"dil_dt", state.dil_dt.data()},
        {"concentration", state.concentration.data()},
        {"source", state.source.data()},
        {"u_x", state.u_x.data()},
        {"u_y", state.u_y.data()},
        {"u_z", state.u_z.data()}
    };
    saver data_saver("suffosion", n_x, n_y, n_z, h);

    //process duration
    const auto duration = n_t * ( (argc > 1) ? std::stoul( argv[1] ) : 1 ); 

	//add wells
    state.wells.insert(std::make_pair(n_x/2, n_y/2));

    prof.tic("time cycle");
    for(auto t = 0; t < duration; ++t)
    {
		solve_pressure(state, ctx, prof);
		
		solve_stress(state, ctx, prof);
		
		solve_chemistry(state);
		
		data_saver.add_step(t * h_t, save_data);
    }
    prof.toc("time cycle");

    std::cout << prof;

    for(auto q = ctx.queue().begin(); q != ctx.queue().end(); ++q)
        q->finish();

    return 0;
}
