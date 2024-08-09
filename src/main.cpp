#include <iostream>
#include "poisson1D.hpp"
#include "solvers.hpp"


#define MAXITER 2000


// @TODO: run with callgrind to ensure no memory leak
int main() {
	Poisson1D problem(
		0.0,
		1.0,
		10,
		[](double x) {return x;},
		{0.0, 0.0}
	);

	// @TODO: add default initialization strategy
	GsSolver solver(problem, InitializationStrategy::Lerp);

	// @BUG: doesn't converge to full precision

	std::cout << "step,rnorm" << std::endl;
	for (int i = 0; i < MAXITER; ++i) {
		std::cout << i << "," << solver.get_residual_norm() << std::endl;

		solver.step();
	}

	return 0;
}
