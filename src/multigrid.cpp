#include "multigrid.hpp"
#include <iostream>


MgSolver::MgSolver(const Poisson1D &problem, const std::vector<MgOp> recipe_, InitializationStrategy strategy) : IterativeSolver(problem, strategy), recipe(recipe_) {
	// very hard memory allocations to understand
}


void MgSolver::step() {

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


bool analyze_cycle_recipe(const std::vector<MgOp> &recipe, int &nlevels) {
	int curr_level = 0;
	int max_depth  = 0;

	for (const auto op : recipe) {
		switch (op) {
			case MgOp::Restrict:
				++curr_level;
				max_depth = std::max(max_depth, curr_level);
				break;

			case MgOp::Prolong:
				--curr_level;

				if (curr_level < 0) {
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

	if (curr_level == 0) {
		return true;
	}
	else {
		std::cerr << "[ERROR]: multigrid cycle doesn't end on the main grid" << std::endl;
		return false;
	}
}
