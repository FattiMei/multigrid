#include "circuit.hpp"
#include "stencil.hpp"
#include "utils.hpp"
#include <stdexcept>


ForcedRC::ForcedRC(
	const double period,
	const double tau_,
	const int n,
	const std::function<double(double)> signal
) :
	Problem(n),
	tau(tau_),
	h(period / (n-1.0)),
	mesh(n)
{
	for (int i = 0; i < n; ++i) {
		mesh[i] = linspace(0.0, period, n, i);
		rhs[i]  = signal(mesh[i]);
	}
}


DiscreteOperator* ForcedRC::get_discrete_operator(const int level) const {
	const double h = std::pow(2.0, level) * this->h;
	int m = n;

	for (int i = 0; i < level; ++i) {
		if ((m-1) % 2 != 0) {
			throw std::invalid_argument("Can't build restricted operator at this level");
		}

		m = 1 + (m-1) / 2;
	}


	// I know it's only a first order formula
	return new ThreePointPeriodicStencil(
		m,
		{-1.0/h + tau, 1.0/h, 0.0}
	);
}


void ForcedRC::set_initial_approximation(double u[], const InitializationStrategy strategy) const {
	(void) strategy;

	for (int i = 0; i < n; ++i) {
		u[i] = 0;
	}
}


const std::vector<double>& ForcedRC::get_mesh() const {
	return mesh;
}
