#include "poisson.hpp"
#include "stencil.hpp"
#include "utils.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>


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

	return new NaiveThreePointStencil(
		m,
		{-1.0 / (h*h), 2.0 / (h*h), -1.0 / (h*h)}
	);
};


DiscreteOperator* PoissonPreciseVariant::get_discrete_operator(const int level) const {
	int m = n;

	for (int i = 0; i < level; ++i) {
		if ((m-1) % 2 != 0) {
			throw std::invalid_argument("Can't build restricted operator at this level");
		}

		m = 1 + (m-1) / 2;
	}

	return new ThreePointStencil(
		m,
		std::pow(2.0, level) * this->h,
		{-1.0, 2.0, -1.0}
	);
};


IsotropicPoisson2D::IsotropicPoisson2D(
	const std::pair<double, double> bottom_left_corner,
	const std::pair<double, double> top_right_corner,
	const int n,
	const std::function<double(double,double)> f,
	const std::function<double(double,double)> boundary
) :
	Problem(n*n),
	rows(n),
	hx((top_right_corner.first - bottom_left_corner.first) / (n-1.0)),
	hy((top_right_corner.second - bottom_left_corner.second) / (n-1.0))
{
	std::vector<double> xx(n);
	std::vector<double> yy(n);

	for (int i = 0; i < n; ++i) {
		xx[i] = linspace(bottom_left_corner.first, top_right_corner.first, n, i);
		yy[i] = linspace(bottom_left_corner.second, top_right_corner.second, n, i);
	}

	for (int i = 0; i < n; ++i) {
		rhs[i] = boundary(xx[i], yy[0]);
	}

	for (int row = 1; row < n; ++row) {
		const int start = n * row;
		const int end   = start + n - 1;

		rhs[start] = boundary(xx[0], yy[row]);

		for (int i = start + 1; i < end; ++i) {
			rhs[i] = f(xx[i - start], yy[row]);
		}

		rhs[end] = boundary(xx[n-1], yy[row]);
	}

	for (int i = 0; i < n; ++i) {
		rhs[n * (n-1) + i] = boundary(xx[i], yy[n-1]);
	}
}


DiscreteOperator* IsotropicPoisson2D::get_discrete_operator(const int level) const {
	const double hx = std::pow(2.0, level) * this->hx;
	const double hy = std::pow(2.0, level) * this->hy;
	int m = rows;

	for (int i = 0; i < level; ++i) {
		if ((m-1) % 2 != 0) {
			throw std::invalid_argument("Can't build restricted operator at this level");
		}

		m = 1 + (m-1) / 2;
	}

	return new FivePointStencil(
		m,
		m,
		{
			/* OVEST  */ -1.0 / (hx*hx),
			/* CENTER */  2.0 / (hx*hx) + 2.0 / (hy*hy),
			/* EST    */ -1.0 / (hx*hx),
			/* NORD   */ -1.0 / (hy*hy),
			/* SUD    */ -1.0 / (hy*hy)
		}
	);
}
