#include "multigrid.hpp"
#include <iostream>
#include <numeric>
#include <cassert>


MgSolver::MgSolver(const Poisson1D &problem, const std::vector<MgOp> recipe_, InitializationStrategy strategy) : IterativeSolver(problem, strategy), recipe(recipe_) {
	if (not analyze_cycle_recipe(recipe, maxlevels)) {
		// throws exception
	}

	if (not compute_grid_sizes(n, maxlevels, grid_size)) {
		// throws exception
	}

	const int total_elements = std::accumulate(grid_size.begin(), grid_size.end(), 0);
	solution_memory = new double[total_elements - n];
	rhs_memory      = new double[total_elements - n];
	residual_memory = new double[total_elements];

	// @DESIGN: this complication is necessary because:
	//   1. It's desirable to do allocations in big chunks (improve heap locality)
	//   2. There is memory already allocated (u comes from BaseSolver, rhs from Poisson1D)
	//
	// grid_solution`, `grid_rhs` and `grid residuals` are maps from level to memory
	// they are tricky to compute but isolate the rest of the program from this complexity
	build_level_to_memory_map();

	for (int level = 0; level < maxlevels+1; ++level) {
		grid_iteration_formula.push_back(problem.get_iteration_formula(level));
		grid_residual_formula.push_back(problem.get_residual_formula(level));
	}
}


MgSolver::~MgSolver() {
	delete[] solution_memory;
	delete[] rhs_memory;
	delete[] residual_memory;
}


void MgSolver::step() {
	int level = 0;

	for (const auto op : recipe) {
		switch (op) {
			case MgOp::Relax: {
				smoother(grid_iteration_formula[level], grid_size[level], grid_rhs[level], grid_solution[level]);
			} break;

			case MgOp::Restrict: {
				// still in doubt about whether to zero the error at the lower level
				const auto formula = grid_residual_formula[level];

				for (int i = 1; i < grid_size[level]; ++i) {
					formula(i, grid_rhs[level], grid_solution[level], grid_residual[level]);
				}

				full_weight_restriction(grid_size[level], grid_residual[level], grid_rhs[level+1]);
				++level;

			} break;


			case MgOp::Prolong: {
				// using residual memory only as alias for the error correction
				double* error_correction = grid_residual[level-1];

				linear_prolongation(grid_size[level], grid_solution[level], error_correction);
				--level;

				for (int i = 1; i < grid_size[level]; ++i) {
					grid_solution[level][i] += error_correction[i];
				}
			} break;


			case MgOp::DirectSolve: {
				assert(0 && "Not implemented");
			} break;


			case MgOp::IterativeSolve: {
				for (int i = 0; i < 300; ++i) {
					smoother(grid_iteration_formula[level], grid_size[level], grid_rhs[level], grid_solution[level]);
				}

			} break;
		}
	}
}


void MgSolver::build_level_to_memory_map() {
	grid_solution.resize(maxlevels+1);
	grid_rhs     .resize(maxlevels+1);
	grid_residual.resize(maxlevels+1);

	// @DESIGN: yet another asymmetry. The rhs at the top level is a `const double *`
	// This pointer is returned by the pointer and it const so that solvers could not overwrite it, changing the problem
	// But now all the rhs pointers are stored in the same vector, to maintain the abstraction I have to break part of it
	//
	// It's either this or make the vector store an std::variant<double*, const double*>
	// This is not really an OOP fault

	grid_solution[0] = u;
	grid_rhs     [0] = rhs;
	grid_residual[0] = residual_memory;

	for (int level = 1; level < maxlevels+1; ++level) {
		grid_residual[level] = grid_residual[level-1] + grid_size[level-1];

		if (level == 1) {
			grid_solution[level] = solution_memory;
			grid_rhs     [level] = rhs_memory;
		}
		else {
			grid_solution[level] = grid_solution[level-1] + grid_size[level-1];
			grid_rhs     [level] = grid_rhs     [level-1] + grid_size[level-1];
		}
	}
}


// all these implementations assume that pointers have sufficient size
void injective_restriction(const int n, const double src[], double dest[]) {
	for (int i = 0, write_index = 0; i < n; i += 2) {
		dest[write_index++] = src[i];
	}
}


void full_weight_restriction(const int n, const double src[], double dest[]) {
	int write_index = 0;
	dest[write_index++] = src[0];

	for (int i = 2; i < n-1; i += 2) {
		dest[write_index++] = 0.25 * src[i-1] + 0.5 * src[i] + 0.25 * src[i+1];
	}

	dest[write_index] = src[n-1];
}


void linear_prolongation(const int m, const double src[], double dest[]) {
	int write_index = 0;

	for (int i = 0; i < m-1; ++i) {
		dest[write_index++] = src[i];
		dest[write_index++] = 0.5 * src[i] + 0.5 * src[i+1];
	}

	dest[write_index] = src[m-1];
}


// the analysis could be more sophisticated but ok
bool analyze_cycle_recipe(const std::vector<MgOp> &recipe, int &nlevels) {
	int level = 0;
	int max_depth  = 0;

	for (const auto op : recipe) {
		switch (op) {
			case MgOp::Restrict:
				++level;
				max_depth = std::max(max_depth, level);
				break;

			case MgOp::Prolong:
				--level;

				if (level < 0) {
					std::cerr << "[ERROR]: multigrid recipe wants to go above finest grid" << std::endl;
					return false;
				}

				break;

			case MgOp::Relax:
			case MgOp::DirectSolve:
			case MgOp::IterativeSolve:
				break;
		}
	}

	nlevels = max_depth;

	if (level == 0) {
		return true;
	}
	else {
		std::cerr << "[ERROR]: multigrid cycle doesn't end on the main grid" << std::endl;
		return false;
	}
}


bool compute_grid_sizes(const int n, const int maxlevels, std::vector<int> &grid_size) {
	int m = n;
	grid_size.push_back(n);

	for (int level = 0; level < maxlevels; ++level) {
		if ((m-1) % 2 == 0) {
			m = 1 + (m-1)/2;

			grid_size.push_back(m);
		}
		else {
			std::cerr << "[ERROR]: can't produce a grid half the size from " << m << " nodes" << std::endl;
			return false;
		}
	}

	return true;
}
