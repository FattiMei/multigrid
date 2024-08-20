#include <iostream>
#include "poisson.hpp"
#include "multigrid.hpp"


int main() {
	constexpr int n       = 2001;
	constexpr int maxiter = 1000;

	IsotropicPoisson1D problem(
		0.0,
		1.0,
		n,
		[](double x){ return x; },
		{0.0, 0.0}
	);

	MgSolver solver(
		&problem,
		MgCycle::V(1),
		InitializationStrategy::Zeros,
		UpdateStrategy::GaussSeidel
	);

	std::cout << "it,mg" << std::endl;

	for (int it = 0; it < maxiter; ++it) {
		std::cout << it << ',';

		std::cout << solver.get_residual_norm();
		solver.step();

		std::cout << std::endl;
	}

	return 0;
}
