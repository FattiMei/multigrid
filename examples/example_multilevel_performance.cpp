#include <iostream>
#include <cmath>
#include <memory>
#include "symbolic.h"
#include "poisson.hpp"
#include "multigrid.hpp"


int main() {
	constexpr int n       = (2 << 10) + 1;
	constexpr int maxiter = 20;

	IsotropicPoisson2D problem(
		{0.0, 0.0},
		{1.0, 1.0},
		n,
		forcing_term_2d,
		solution_2d
	);

	const auto smoother     = UpdateStrategy::GaussSeidel;
	const auto restriction  = full_weight_restriction_2d;
	const auto prolongation = linear_prolongation_2d;

	// when problems are so big the direct solver stalls, so no 2 level
	const std::vector<std::pair<std::string,std::shared_ptr<IterativeSolver>>> solvers{
		// {"2-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(2, 3, false), InitializationStrategy::Zeros, smoother, restriction, prolongation)},
		// {"3-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(3, 3), InitializationStrategy::Zeros, smoother, restriction, prolongation)},
		{"V(5)-GS"   , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, smoother, restriction, prolongation)},
		{"V(7) GS"   , std::make_shared<MgSolver>(&problem, MgCycle::V(7, 3), InitializationStrategy::Zeros, smoother, restriction, prolongation)}
	};

	std::cout << "it";
	for (auto [label, _] : solvers) {
		std::cout << "," << label;
	}
	std::cout << std::endl;


	for (int i = 0; i < maxiter; ++i) {
		std::cout << i;

		for (auto& [_, solver] : solvers) {
			std::cout << "," << solver->get_residual_norm();
			solver->step();
		}

		std::cout << std::endl;
	}


	return 0;
}
