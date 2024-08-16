#include <iostream>
#include "poisson1D.hpp"
#include "solvers.hpp"
#include "multigrid.hpp"


#define MAXITER 2000


// @TODO: run with callgrind to ensure no memory leak
int main() {
	IsotropicPoisson1D problem2(
		0.0,
		1.0,
		10,
		[](double x) {return x;},
		{0.0, 0.0}
	);

	return 0;
}
