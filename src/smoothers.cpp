#include "smoothers.hpp"


Smoother::Smoother(const Poisson1D &problem) : iteration_formula(problem.get_iteration_formula()) {}


Jacobi::Jacobi(const Poisson1D &problem) : Smoother(problem) {
	// allocate local working memory
}


Jacobi::~Jacobi() {
	// deallocate local working memory
}


double* Jacobi::smooth(const int n, const double b[], double u[]) {
	(void) n;
	(void) b;

	return u;
}


double* GS::smooth(const int n, const double b[], double u[]) {
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
