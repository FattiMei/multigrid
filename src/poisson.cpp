#include "poisson.hpp"
#include "utils.hpp"
#include <cmath>
#include <stdexcept>


IsotropicPoisson1D::IsotropicPoisson1D(
	const double inf,
	const double sup,
	const int    n,
	const std::function<double(double)> f,
	const std::pair<double,double> boundary
) :
	Problem(n),
	h((sup - inf) / (n - 1.0))
{
	rhs[0] = boundary.first;

	for (int i = 1; i < n-1; ++i) {
		rhs[i] = f(linspace(inf, sup, n, i));
	}

	rhs[n-1] = boundary.second;
}


DiscreteOperator* IsotropicPoisson1D::get_discrete_operator(const int level) const {
	const double h = std::pow(2.0, level) * this->h;
	int m = n;

	for (int i = 0; i < level; ++i) {
		if ((m-1) % 2 != 0) {
			throw std::invalid_argument("Can't build restricted operator at this level");
		}

		m = 1 + (m-1) / 2;
	}

	return new ThreePointStencil(
		m,
		{-1.0 / (h*h), 2.0 / (h*h), -1.0 / (h*h)}
	);
};
