#include <iostream>
#include <memory>
#include "poisson.hpp"
#include "solvers.hpp"


int main() {
	constexpr int n       = 100;
	constexpr int maxiter = 1000;

	IsotropicPoisson1D problem(
		0.0,
		1.0,
		n,
		[](double x){ return x; },
		{0.0, 0.0}
	);


	const std::vector<std::pair<std::string, std::shared_ptr<IterativeSolver>>> solvers{
		{"jacobi"   , std::make_shared<SmootherSolver>(&problem, InitializationStrategy::Zeros, UpdateStrategy::Jacobi)},
		{"gseidel"  , std::make_shared<SmootherSolver>(&problem, InitializationStrategy::Zeros, UpdateStrategy::GaussSeidel)},
		{"red-black", std::make_shared<SmootherSolver>(&problem, InitializationStrategy::Zeros, UpdateStrategy::RedBlack)}
	};


	std::cout << "n";
	for (auto [label, _] : solvers) std::cout << ',' << label;
	std::cout << std::endl;


	for (int it = 0; it <= maxiter; ++it) {
		std::cout << it;

		for (auto [_, solver] : solvers) {
			std::cout << ',' << solver->get_residual_norm();
			solver->step();
		}

		std::cout << std::endl;
	}


	return 0;
}
