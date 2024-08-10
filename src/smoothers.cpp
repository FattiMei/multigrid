#include "smoothers.hpp"


using namespace Smoother;


BaseSmoother::BaseSmoother(const Poisson1D &problem) : iteration_formula(problem.get_iteration_formula()) {}


Jacobi::Jacobi(const Poisson1D &problem) : BaseSmoother(problem) {
	local = new double[problem.get_problem_size()];
}


Jacobi::~Jacobi() {
	delete[] local;
}


double* Jacobi::smooth(const int n, const double b[], double u[]) {
	// Jersey style, does two sweeps
	for (int i = 1; i < n-1; ++i) {
		iteration_formula(i, b, u, local);
	}

	for (int i = 1; i < n-1; ++i) {
		iteration_formula(i, b, local, u);
	}

	return u;
}


double* GSeidel::smooth(const int n, const double b[], double u[]) {
	for (int i = 1; i < n-1; ++i) {
		iteration_formula(i, b, u, u);
	}

	return u;
}


double* RedBlack::smooth(const int n, const double b[], double u[]) {
	for (int i = 1; i < n-1; i += 2) {
		iteration_formula(i, b, u, u);
	}

	for (int i = 2; i < n-1; i += 2) {
		iteration_formula(i, b, u, u);
	}

	return u;
}
