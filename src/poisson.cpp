#include "poisson.hpp"
#include "stencil.hpp"
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

	return new NaiveThreePointStencil(
		m,
		{-1.0 / (h*h), 2.0 / (h*h), -1.0 / (h*h)}
	);
}


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
}


void IsotropicPoisson1D::set_initial_approximation(double u[], const InitializationStrategy strategy) const {
	u[0] = rhs[0];
	u[n-1] = rhs[n-1];

	switch (strategy) {
		case InitializationStrategy::Zeros: {
			for (int i = 1; i < n-1; ++i) {
				u[i] = 0.0;
			}
		} break;

		case InitializationStrategy::Lerp: {
			for (int i = 1; i < n-1; ++i) {
				u[i] = linspace(rhs[0], rhs[n-1], n, i);
			}
		} break;
	}
}


IsotropicPoisson2D::IsotropicPoisson2D(
	const std::pair<double, double> bottom_left_corner,
	const std::pair<double, double> top_right_corner,
	const int n,
	const std::function<double(double,double)> f,
	const std::function<double(double,double)> boundary
) :
	Problem(n*n),
	rows(n),
	cols(n),
	hx((top_right_corner.first  - bottom_left_corner.first) / (n-1.0)),
	hy((top_right_corner.second - bottom_left_corner.second) / (n-1.0)),
	xx(n),
	yy(n)
{
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


AnisotropicPoisson2D::AnisotropicPoisson2D(
	const std::pair<double, double> bottom_left_corner,
	const std::pair<double, double> top_right_corner,
	const int n,
	const std::function<double(double,double)> f,
	const std::function<double(double,double)> boundary,
	const std::function<double(double,double)> c
) :
	IsotropicPoisson2D(bottom_left_corner, top_right_corner, n, f, boundary)
{
	// @DESIGN: I'm aware of the repetition in mesh calculation and rhs access, but it's partially justified by the one time use
	// To solve this performance problem I would have to make this the default class and the IsotropicPoisson2D would be just constructing this class
	// with c(x,y) == 1

	for (int row = 1; row < n; ++row) {
		const int start = n * row;
		const int end   = start + n - 1;

		for (int i = start + 1; i < end; ++i) {
			rhs[i] /= c(xx[i - start], yy[row]);
		}
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


void IsotropicPoisson2D::set_initial_approximation(double u[], const InitializationStrategy strategy) const {
	// ignoring the lerp strategy, future work might expand on this
	(void) strategy;

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			if (i == 0 or i == (rows-1) or j == 0 or j == (cols-1)) {
				u[i*cols+j] = rhs[i*cols+j];
			}
			else {
				u[i*cols+j] = 0.0;
			}
		}
	}
}
