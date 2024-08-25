#ifndef __POISSON_HPP__
#define __POISSON_HPP__


#include "problem.hpp"


class IsotropicPoisson1D : public Problem {
	public:
		IsotropicPoisson1D(
			const double inf,
			const double sup,
			const int    n,
			const std::function<double(double)> f,
			const std::pair<double,double> boundary
		);
		~IsotropicPoisson1D() = default;

		double get_step() const override { return h; };
		DiscreteOperator* get_discrete_operator(const int level = 0) const override;


	private:
		const double h;
};


class PoissonPreciseVariant : public Problem {
	public:
		PoissonPreciseVariant(
			const double inf,
			const double sup,
			const int    n,
			const std::function<double(double)> f,
			const std::pair<double,double> boundary
		);
		~PoissonPreciseVariant() = default;

		double get_step() const override { return h; };
		DiscreteOperator* get_discrete_operator(const int level = 0) const override;


	private:
		const double h;
};


class IsotropicPoisson2D : public Problem {
	public:
		IsotropicPoisson2D(
			const std::pair<double, double> bottom_left_corner,
			const std::pair<double, double> top_right_corner,
			const int n,
			const std::function<double(double,double)> f,
			const std::function<double(double,double)> boundary
		);
		~IsotropicPoisson2D() = default;

		double get_step() const override { return hx; };
		DiscreteOperator* get_discrete_operator(const int level = 0) const override;


	private:
		const int rows;
		const double hx;
		const double hy;
};


#endif
