#include "circuit.hpp"
#include "stencil.hpp"
#include "utils.hpp"
#include <stdexcept>


ForcedRC::ForcedRC(
	const double period,
	const double RC,
	const int n,
	const std::function<double(double)> signal
) :
	Problem(n),
	tau(RC),
	h(period / n),
	mesh(n)
{
	for (int i = 0; i < n; ++i) {
		mesh[i] = linspace(0.0, period, n+1, i);
		rhs[i]  = signal(mesh[i]);
	}
}


DiscreteOperator* ForcedRC::get_discrete_operator(const int level) const {
	if (level != 0) {
		throw std::invalid_argument("Admits only fine grid operator");
	}

	return new ThreePointPeriodicStencil(
		n,
		{-1.0/h + tau, 1.0/h, 0.0}
	);
}


const std::vector<double>& ForcedRC::get_mesh() const {
	return mesh;
}
