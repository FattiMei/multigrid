#include "smoothers.hpp"


using namespace Smoother;


BaseSmoother::BaseSmoother(const int n_, const Update iteration_formula_) : iteration_formula(iteration_formula_), n(n_) {}


Jacobi::Jacobi(const int n, const Update iteration_formula) : BaseSmoother(n, iteration_formula) {
	local = new double[n];
}


Jacobi::~Jacobi() {
	delete[] local;
}


void Jacobi::smooth(const double b[], double u[]) {
	// Jersey style, does two sweeps
	for (int i = 1; i < n-1; ++i) {
		iteration_formula(i, b, u, local);
	}

	for (int i = 1; i < n-1; ++i) {
		iteration_formula(i, b, local, u);
	}
}


void GSeidel::smooth(const double b[], double u[]) {
	for (int i = 1; i < n-1; ++i) {
		iteration_formula(i, b, u, u);
	}
}


void RedBlack::smooth(const double b[], double u[]) {
	for (int i = 1; i < n-1; i += 2) {
		iteration_formula(i, b, u, u);
	}

	for (int i = 2; i < n-1; i += 2) {
		iteration_formula(i, b, u, u);
	}
}
