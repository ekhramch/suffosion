#pragma once
#include <vector>
#include <vexcl/vexcl.hpp>
#include <amgcl/profiler.hpp>
#include "utils.hpp"

//solve pressure equation
void solve_pressure(Unknowns &state, vex::Context &ctx, amgcl::profiler &prof);

//build matrix for pressure equation
void assembly_matrix(
						std::vector<int> &col, 
						std::vector<double> &val, 
						std::vector<int> &ptr, 
						std::vector<double> &rhs, 
						Unknowns &state
					);