#include <iostream>
#include <cmath>
#include <memory>
#include "poisson.hpp"
#include "multigrid.hpp"


int main() {
	constexpr int n       = (2 << 8) + 1;
	constexpr int maxiter = 10;

	IsotropicPoisson2D problem(
		{0.0, 0.0},
		{1.0, 1.0},
		n,
		[](double x, double y){ return x + y; },
		[](double x, double y){ return 0.0; }
	);

	const auto smoother     = UpdateStrategy::GaussSeidel;
	const auto restriction  = injection_restriction_2d;
	const auto prolongation = linear_prolongation_2d;

	const std::vector<std::pair<std::string,std::shared_ptr<IterativeSolver>>> solvers{
		// {"2-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(2, 3, false), InitializationStrategy::Zeros, smoother, restriction, prolongation)},
		{"3-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(3, 3), InitializationStrategy::Zeros, smoother, restriction, prolongation)},
		{"5-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(5, 3), InitializationStrategy::Zeros, smoother, restriction, prolongation)},
		{"7-level"   , std::make_shared<MgSolver>(&problem, MgCycle::V(7, 3), InitializationStrategy::Zeros, smoother, restriction, prolongation)}
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
