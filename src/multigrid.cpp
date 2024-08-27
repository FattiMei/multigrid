#include "multigrid.hpp"
#include <iostream>
#include <numeric>


std::ostream& operator<<(std::ostream& os, const MgOp op) {
	switch (op) {
		case MgOp::Relax         : os << "Relax"         ; break;
		case MgOp::Restrict      : os << "Restrict"      ; break;
		case MgOp::Prolong       : os << "Prolong"       ; break;
		case MgOp::DirectSolve   : os << "DirectSolve"   ; break;
		case MgOp::IterativeSolve: os << "IterativeSolve"; break;
	}

	return os;
}


MgSolver::MgSolver(
	const Problem*			problem,
	const std::vector<MgOp>		cycle_spec,
	const InitializationStrategy	strategy,
	const UpdateStrategy		smoother_of_choice,
	const RestrictionOperator	restrictor,
	const ProlongationOperator	prolonger
) :
	IterativeSolver(problem, strategy),
	n(problem->get_size()),
	recipe(cycle_spec),
	smoother(smoother_of_choice),
	maxlevels(analyze_cycle_recipe(recipe)),
	restrict(restrictor),
	prolong(prolonger),
	prova(compute_grid_sizes(problem, maxlevels)),
	grid_size(compute_grid_sizes(n, maxlevels))
{
	const int total_elements = std::accumulate(grid_size.begin(), grid_size.end(), 0);
	solution_memory = new double[total_elements - n];
	rhs_memory      = new double[total_elements - n];
	residual_memory = new double[total_elements];


	// @DESIGN: this complication is necessary because:
	//   1. It's desirable to do allocations in big chunks (improve heap locality)
	//   2. There is memory already allocated (u comes from BaseSolver, rhs from Poisson1D)
	//
	// `grid_solution`, `grid_rhs` and `grid residuals` are maps from level to memory
	// they are tricky to compute but isolate the rest of the program from this complexity


	build_level_to_memory_map();

	grid_operator.push_back(op);
	for (int level = 1; level <= maxlevels; ++level) {
		grid_operator.push_back(problem->get_discrete_operator(level));
	}

	direct_solver.compute(grid_operator.back()->get_sparse_repr());
}


MgSolver::~MgSolver() {
	delete[] solution_memory;
	delete[] rhs_memory;
	delete[] residual_memory;

	for (int level = 1; level <= maxlevels; ++level) {
		delete grid_operator[level];
	}
}


void MgSolver::step() {
	int level = 0;

	for (const auto op : recipe) {
		switch (op) {
			case MgOp::Relax: {
				grid_operator[level]->relax(grid_rhs[level].get_const_ptr(), grid_solution[level], smoother);
			} break;

			case MgOp::Restrict: {
				grid_operator[level]->compute_residual(grid_rhs[level].get_const_ptr(), grid_solution[level], grid_residual[level]);

				// zeroing the error is important!
				for (int i = 0; i < grid_size[level+1]; ++i) grid_solution[level+1][i] = 0.0;

				restrict(grid_size[level], grid_residual[level], grid_rhs[level+1].get_mutable_ptr());
				++level;
			} break;

			case MgOp::Prolong: {
				// using residual memory only as alias for the error correction
				double* error_correction = grid_residual[level-1];

				prolong(grid_size[level], grid_solution[level], error_correction);
				--level;

				for (int i = 0; i < grid_size[level]; ++i) {
					grid_solution[level][i] += error_correction[i];
				}
			} break;

			case MgOp::DirectSolve: {
				Eigen::Map<const Eigen::VectorXd> b(grid_rhs[level].get_const_ptr(), grid_size[level]);
				Eigen::Map<Eigen::VectorXd> tmp(grid_solution[level], grid_size[level]);

				tmp = direct_solver.solve(b);
			} break;

			case MgOp::IterativeSolve: {
				for (int i = 0; i < 300; ++i) {
					grid_operator[level]->relax(grid_rhs[level].get_const_ptr(), grid_solution[level], smoother);
				}
			} break;
		}
	}
}


void MgSolver::build_level_to_memory_map() {
	grid_solution.resize(maxlevels+1);
	grid_residual.resize(maxlevels+1);
	grid_rhs     .resize(maxlevels+1);

	grid_solution[0] = u.data();
	grid_rhs     [0] = PointerVariant<double>(rhs);
	grid_residual[0] = residual_memory;

	for (int level = 1; level < maxlevels+1; ++level) {
		grid_residual[level] = grid_residual[level-1] + grid_size[level-1];

		if (level == 1) {
			grid_solution[level] = solution_memory;
			grid_rhs     [level] = PointerVariant<double>(rhs_memory);
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
int analyze_cycle_recipe(const std::vector<MgOp> &recipe) {
	int level    = 0;
	int maxdepth = 0;

	for (const auto op : recipe) {
		switch (op) {
			case MgOp::Restrict:
				++level;
				maxdepth = std::max(maxdepth, level);
				break;

			case MgOp::Prolong:
				--level;

				if (level < 0) {
					throw std::logic_error("cycle spec doesn't end on finest grid");
				}

				break;

			case MgOp::Relax:
			case MgOp::DirectSolve:
			case MgOp::IterativeSolve:
				break;
		}
	}

	return maxdepth;
}


std::vector<int> compute_grid_sizes(const int n, const int maxdepth) {
	std::vector<int> grid_size(maxdepth+1);
	int m = n;

	grid_size[0] = n;

	for (int level = 1; level <= maxdepth; ++level) {
		if ((m-1) % 2 == 0) {
			m = 1 + (m-1)/2;

			grid_size[level] = m;
		}
		else {
			throw std::logic_error("Can't produce subgrid");
		}
	}

	return grid_size;
}


std::vector<std::pair<int,int>> compute_grid_sizes(const Problem* problem, const int maxdepth) {
	std::vector<std::pair<int,int>> grid_size(maxdepth+1);
	grid_size[0] = {problem->get_dimension(0), problem->get_dimension(1)};

	for (int level = 1; level <= maxdepth; ++level) {
		int m = grid_size[0].first;

		if ((m-1) % 2 == 0) {
			m = 1 + (m-1)/2;

			grid_size[level].first = m;
		}
		else {
			throw std::logic_error("Can't produce subgrid");
		}
	}

	if (grid_size[0].second != 0) {
		int m = grid_size[0].second;

		for (int level = 1; level <= maxdepth; ++level) {
			if ((m-1) % 2 == 0) {
				m = 1 + (m-1)/2;

				grid_size[level].second = m;
			}
			else {
				throw std::logic_error("Can't produce subgrid");
			}
		}
	}

	/* DEBUG
	for (auto& [a, b] : grid_size) {
		std::cout << a << ' ' << b << std::endl;
	}
	*/

	return grid_size;
}


// I could in principle precompute the number of MgOps given the cycle spec to optimize the allocations...
std::vector<MgOp> MgCycle::V(const int levels_one_indexed, const int smoothing_steps, const bool solve) {
	const int levels = levels_one_indexed - 1;
	std::vector<MgOp> result;

	for (int i = 0; i < levels; ++i) {
		for (int step = 0; step < smoothing_steps; ++step) {
			result.push_back(MgOp::Relax);
		}

		result.push_back(MgOp::Restrict);
	}

	if (solve) {
		result.push_back(MgOp::DirectSolve);
	}
	else {
		result.push_back(MgOp::Relax);
	}

	for (int i = 0; i < levels; ++i) {
		result.push_back(MgOp::Prolong);

		for (int step = 0; step < smoothing_steps; ++step) {
			result.push_back(MgOp::Relax);
		}
	}

	return result;
}


std::vector<MgOp> MgCycle::F(const int levels_one_indexed, const int smoothing_steps, const bool solve) {
	const int levels   = levels_one_indexed - 1;
	const auto MgSolve = solve ? MgOp::DirectSolve : MgOp::Relax;
	std::vector<MgOp> result;

	for (int i = 0; i < levels; ++i) {
		for (int step = 0; step < smoothing_steps; ++step) {
			result.push_back(MgOp::Relax);
		}

		result.push_back(MgOp::Restrict);
	}

	for (int i = 0; i < levels; ++i) {
		for (int j = 0; j < i; ++j) {
			result.push_back(MgOp::Prolong);

			for (int step = 0; step < smoothing_steps; ++step) {
				result.push_back(MgOp::Relax);
			}
		}

		for (int j = 0; j < i; ++j) {
			result.push_back(MgOp::Restrict);

			if (j < (i-1)) {
				for (int step = 0; step < smoothing_steps; ++step) {
					result.push_back(MgOp::Relax);
				}
			}
		}

		result.push_back(MgSolve);
	}
	
	for (int i = 0; i < levels; ++i) {
		result.push_back(MgOp::Prolong);

		for (int step = 0; step < smoothing_steps; ++step) {
			result.push_back(MgOp::Relax);
		}
	}

	return result;
}


void wcycle_helper(std::vector<MgOp>& recipe, const int level, const int maxlevels, const int smoothing_steps, const MgOp MgSolve) {
	if (level == maxlevels) {
		for (int i = 0; i < smoothing_steps; ++i) {
			recipe.push_back(MgOp::Relax);
		}

		recipe.push_back(MgOp::Restrict);
		recipe.push_back(MgSolve);
		recipe.push_back(MgOp::Prolong);
	}
	else {
		for (int i = 0; i < smoothing_steps; ++i) {
			recipe.push_back(MgOp::Relax);
		}
		recipe.push_back(MgOp::Restrict);

		wcycle_helper(recipe, level+1, maxlevels, smoothing_steps, MgSolve);
		wcycle_helper(recipe, level+1, maxlevels, smoothing_steps, MgSolve);

		for (int i = 0; i < smoothing_steps; ++i) {
			recipe.push_back(MgOp::Relax);
		}
		recipe.push_back(MgOp::Prolong);
	}
}


std::vector<MgOp> MgCycle::W(const int levels_one_indexed, const int smoothing_steps, const bool solve) {
	const int levels   = levels_one_indexed - 1;
	const auto MgSolve = solve ? MgOp::DirectSolve : MgOp::Relax;
	std::vector<MgOp> result;

	// recursive approach to the problem
	wcycle_helper(result, 1, levels, smoothing_steps, MgSolve);
	result.push_back(MgOp::Relax);

	return result;
}
