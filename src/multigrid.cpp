#include "multigrid.hpp"
#include <iostream>
#include <numeric>
#include <algorithm>


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
	grid_dim(compute_grid_dim(problem, maxlevels)),
	grid_size(compute_grid_size(grid_dim))
{
	// @DESIGN: zero initialize the memory so that projection and restriction operators don't need to care about the boundary nodes
	const int total_elements = std::accumulate(grid_size.begin(), grid_size.end(), 0);
	solution_memory = new double[total_elements - n]{};
	rhs_memory      = new double[total_elements - n]{};
	residual_memory = new double[total_elements]{};


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

	// for some reason I couldn't do it in the function `analyze_cycle_recipe` in the constructor
	if (std::find(recipe.begin(), recipe.end(), MgOp::DirectSolve) != std::end(recipe)) {
		Eigen::SparseMatrix<double> A = grid_operator.back()->get_sparse_repr();
		direct_solver.analyzePattern(A);
		direct_solver.factorize(A);

		if (direct_solver.info() != Eigen::Success) {
			throw std::runtime_error("Eigen has failed to factorize the matrix");
		}

		// direct_solver.compute(grid_operator.back()->get_sparse_repr());
	}
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

				restrict(grid_dim[level], grid_residual[level], grid_rhs[level+1].get_mutable_ptr());
				++level;
			} break;

			case MgOp::Prolong: {
				// using residual memory only as alias for the error correction
				double* error_correction = grid_residual[level-1];

				prolong(grid_dim[level-1], grid_solution[level], error_correction);
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
// using index arithmetic instead of write_index to avoid loop carried dependencies
void injection_restriction_1d(const std::pair<int,int> dim, const double src[], double dest[]) {
	for (int i = 0; i < dim.first; i += 2) {
		dest[i / 2] = src[i];
	}
}


void full_weight_restriction_1d(const std::pair<int,int> dim, const double src[], double dest[]) {
	for (int i = 2; i < dim.first-1; i += 2) {
		dest[i / 2] = 0.25 * src[i-1] + 0.5 * src[i] + 0.25 * src[i+1];
	}
}


// dim is the TARGET grid dimension
void linear_prolongation_1d(const std::pair<int,int> dim, const double src[], double dest[]) {
	for (int i = 0; i < dim.first; ++i) {
		if (i % 2 == 0) {
			dest[i] = src[i / 2];
		}
		else {
			dest[i] = 0.5 * src[i / 2] + 0.5 * src[i / 2 + 1];
		}
	}
}


void injection_restriction_2d(const std::pair<int,int> dim, const double src[], double dest[]) {
	const int target_cols = 1 + (dim.second - 1) / 2;

	for (int i = 0; i < dim.first; i += 2) {
		for (int j = 0; j < dim.second; j += 2) {
			dest[(i / 2) * target_cols + (j / 2)] = src[i * dim.second + j];
		}
	}
}


void full_weight_restriction_2d(const std::pair<int,int> dim, const double src[], double dest[]) {
	const int target_cols = 1 + (dim.second - 1) / 2;

	for (int i = 2; i < dim.first-1; i += 2) {
		for (int j = 2; j < dim.second; j += 2) {
			const int linear_index = i * dim.second + j;

			/*
			dest[(i / 2) * target_cols + (j / 2)] =
				  0.5   * src[linear_index]
				+ 0.125 * src[linear_index - 1]
				+ 0.125 * src[linear_index + 1]
				+ 0.125 * src[linear_index + dim.second]
				+ 0.125 * src[linear_index - dim.second];
			*/

			dest[(i / 2) * target_cols + (j / 2)] =
				  0.25  * src[linear_index]
				+ 0.125 * src[linear_index - 1]
				+ 0.125 * src[linear_index + 1]
				+ 0.125 * src[linear_index + dim.second]
				+ 0.125 * src[linear_index - dim.second]
				+ 0.0625 * src[linear_index + dim.second + 1]
				+ 0.0625 * src[linear_index + dim.second - 1]
				+ 0.0625 * src[linear_index - dim.second + 1]
				+ 0.0625 * src[linear_index - dim.second - 1];
		}
	}
}


void linear_prolongation_2d(const std::pair<int,int> dim, const double src[], double dest[]) {
	const int source_cols = 1 + (dim.second - 1) / 2;

	for (int i = 1; i < dim.first-1; ++i) {
		for (int j = 1; j < dim.second-1; ++j) {
			const int linear_index = (i / 2) * source_cols + (j / 2);

			if (i % 2 == 1 and j % 2 == 1) {
				dest[i * dim.second + j] =
					  0.25 * src[linear_index]
					+ 0.25 * src[linear_index + 1]
					+ 0.25 * src[linear_index + source_cols]
					+ 0.25 * src[linear_index + source_cols + 1];
			}
			else if (i % 2 == 1) {
				dest[i * dim.second + j] =
					  0.5 * src[linear_index]
					+ 0.5 * src[linear_index + source_cols];
			}
			else if (j % 2 == 1) {
				dest[i * dim.second + j] =
					  0.5 * src[linear_index]
					+ 0.5 * src[linear_index + 1];
			}
			else {
				dest[i * dim.second + j] = src[linear_index];
			}
		}
	}
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

			case MgOp::DirectSolve:
			case MgOp::Relax:
			case MgOp::IterativeSolve:
				break;
		}
	}

	return maxdepth;
}


std::vector<std::pair<int,int>> compute_grid_dim(const Problem* problem, const int maxdepth) {
	std::vector<std::pair<int,int>> grid_size(maxdepth+1);
	grid_size[0] = {problem->get_dimension(0), problem->get_dimension(1)};

	int m = grid_size[0].first;
	for (int level = 1; level <= maxdepth; ++level) {

		if ((m-1) % 2 == 0) {
			m = 1 + (m-1)/2;

			grid_size[level].first = m;
		}
		else {
			throw std::logic_error("Can't produce subgrid");
		}
	}

	m = grid_size[0].second;
	if (m != 0) {
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


std::vector<int> compute_grid_size(const std::vector<std::pair<int,int>>& sizes) {
	std::vector<int> result(sizes.size());

	for (size_t i = 0; i < result.size(); ++i) {
		result[i] = sizes[i].first;

		if (sizes[i].second != 0) {
			result[i] *= sizes[i].second;
		}
	}

	return result;
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
