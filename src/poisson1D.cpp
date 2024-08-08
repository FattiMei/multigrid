#include "poisson1D.hpp"
#include "utils.hpp"
#include <cmath>


Poisson1D::Poisson1D(double inf, double sup, int n_, std::function<double(double)> f, std::pair<double,double> boundary) : n(n_), h((sup - inf) / (n-1.0)) {
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


void Poisson1D::set_initial_approximation(double *u, InitializationStrategy strategy) const {
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


int Poisson1D::get_problem_size() const {
	return n;
}


const double* Poisson1D::get_rhs() const {
	return static_cast<const double *>(b);
}


const Update Poisson1D::get_iteration_formula() const {
	return [this](int i, const double b[], const double src[], double dest[]) {
		const double h = this->h;

		dest[i] = (h*h*b[i] + src[i-1] + src[i+1]) / 2.0;
	};
}


double Poisson1D::get_residual_norm(const double u[]) const {
	double norm = 0.0;

	for (int i = 1; i < n-1; ++i) {
		norm += b[i] - (2.0*u[i] - u[i-1] - u[i+1]) / (h*h);
	}

	return std::sqrt(norm);
}
