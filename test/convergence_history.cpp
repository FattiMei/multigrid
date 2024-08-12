#include <iostream>
#include <vector>
#include "poisson1D.hpp"
#include "smoothers.hpp"
#include "solvers.hpp"
#include "multigrid.hpp"


// @TODO(minor detail) make the problem size a command line parameter
int main() {
	constexpr int n       = 1001;
	constexpr int maxiter = 2000;

	const Poisson1D problem(
		0.0,
		1.0,
		n,
		[](double x) {return x;},
		{0.0, 0.0}
	);


	const std::vector<MgOp> recipe {
		MgOp::Relax,
		MgOp::Restrict,
		MgOp::Relax,
		MgOp::Prolong,
		MgOp::Relax
	};


	const std::vector<MgOp> only_smooth {
		MgOp::Relax,
		MgOp::Relax
	};


	// nice declarative style, a little too cumbersome IMO
	const std::vector<std::pair<std::string,IterativeSolver*>> solvers{
		{"gseidel"     , new SmootherSolver<Smoother::GSeidel> (problem, InitializationStrategy::Zeros)},
		{"gseidel^2"   , new MgSolver(problem, only_smooth, InitializationStrategy::Zeros)},
		{"two-level MG", new MgSolver(problem, recipe, InitializationStrategy::Zeros)},
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
	// @TODO: use callgrind to assess no memory leak
	for (auto [_, solver] : solvers) {
		delete solver;
	}


	return 0;
}
