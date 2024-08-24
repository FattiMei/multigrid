#include <iostream>
#include <memory>
#include "poisson.hpp"
#include "multigrid.hpp"


int main() {
	constexpr int n       = (2 << 15) + 1;
	constexpr int maxiter = 100;

	const IsotropicPoisson1D problem(
		0.0,
		1.0,
		n,
		[](double x){ return x; },
		{0.0, 0.0}
	);

	// TODO: use std::unique_ptr
	const std::vector<std::pair<std::string,std::shared_ptr<IterativeSolver>>> solvers{
		{"V-cycle"   , std::make_shared<MgSolver>(&problem, MgCycle::V(4), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction, linear_prolongation)},
		{"F-cycle"   , std::make_shared<MgSolver>(&problem, MgCycle::F(4), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction, linear_prolongation)},
		{"W-cycle"   , std::make_shared<MgSolver>(&problem, MgCycle::W(4), InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel, full_weight_restriction, linear_prolongation)},
	};

	std::cout << "it";
	for (auto [label, _] : solvers) {
		std::cout << "," << label;
	}
	std::cout << std::endl;


	for (int i = 0; i < maxiter; ++i) {
		std::cout << i;

		for (auto [_, solver] : solvers) {
			std::cout << "," << solver->get_residual_norm();
			solver->step();
		}

		std::cout << std::endl;
	}


	return 0;
}
