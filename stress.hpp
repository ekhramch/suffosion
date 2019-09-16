#pragma once
#include <vector>
#include <set>
#include <utilities>
#include <vexcl/vexcl.hpp>
#include <amgcl/profiler.hpp>
#include "utils.hpp"

//solve Lame system
void solve_stress(Unknowns &state, vex::Context &ctx, amgcl::profiler &prof);

//build matrix for Lame system
void assembly_matrix(
						std::vector<int> &col, 
						std::vector<double> &val, 
						std::vector<int> &ptr,
						std::vector<double> &rhs,												
						std::vector<double> &pressure, 
						std::set<std::pair<int, int>> wells
					);

//calculate trace of strain tensor
void get_dilatation(Unknowns &state);
