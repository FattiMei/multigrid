#include <iostream>
#include <vector>
#include "poisson1D.hpp"
#include "smoothers.hpp"
#include "solvers.hpp"


// @TODO(minor detail) make the problem size a command line parameter
int main() {
	constexpr int n       = 1000;
	constexpr int maxiter = 2000;

	const Poisson1D problem(
		0.0,
		1.0,
		n,
		[](double x) {return x;},
		{0.0, 0.0}
	);


	// nice declarative style, a little too cumbersome IMO
	const std::vector<std::pair<std::string,IterativeSolver*>> solvers{
		{"gseidel", new SmootherSolver<Smoother::GSeidel>(problem, InitializationStrategy::Zeros)}
	};


	// csv header
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


	// good practice to free memory
	for (auto [_, solver] : solvers) {
		delete solver;
	}


	return 0;
}
