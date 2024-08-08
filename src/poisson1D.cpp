#include "poisson1D.hpp"
#include "utils.hpp"


Poisson1D::Poisson1D(double inf, double sup, int n_, std::function<double(double)> f, std::pair<double,double> boundary) : n(n_) {
	b = new double[n];
	b[0] = boundary.first;

	for (int i = 1; i < n-1; ++i) {
		b[i] = f(linspace(inf, sup, n, i));
	}

	b[n-1] = boundary.second;
}


Poisson1D::~Poisson1D() {
	delete[] b;
}


void Poisson1D::set_initial_approximation(double *u, InitializationStrategy strategy) {
	u[0] = b[0];

	switch(strategy) {
		case InitializationStrategy::Zeros:
			for (int i = 1; i < n-1; ++i) {
				u[i] = 0.0;
			}

			break;

		case InitializationStrategy::Lerp:
			for (int i = 1; i < n-1; ++i) {
				u[i] = linspace(b[0], b[n-1], n, i);
			}

			break;
	}

	u[n-1] = b[n-1];
}


int Poisson1D::get_problem_size() {
	return n;
}


const double* Poisson1D::get_rhs() {
	return static_cast<const double *>(b);
}
