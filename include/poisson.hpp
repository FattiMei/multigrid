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
		int get_dimension(const int dim) const override {
			return dim == 0 ? n : 0;
		};
		DiscreteOperator* get_discrete_operator(const int level = 0) const override;

		void set_initial_approximation(double u[], const InitializationStrategy strategy) const override;


	protected:
		const double h;
};


class PoissonPreciseVariant : public IsotropicPoisson1D {
	public:
		PoissonPreciseVariant(
			const double inf,
			const double sup,
			const int    n,
			const std::function<double(double)> f,
			const std::pair<double,double> boundary
		) : IsotropicPoisson1D(inf, sup, n, f, boundary) {};

		DiscreteOperator* get_discrete_operator(const int level = 0) const override;
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
		int get_dimension(const int dim) const override {
			if (dim == 0) {
				return rows;
			}

			if (dim == 1) {
				return cols;
			}

			return 0;
		}
		DiscreteOperator* get_discrete_operator(const int level = 0) const override;

		void set_initial_approximation(double u[], const InitializationStrategy strategy) const override;


	protected:
		const int rows;
		const int cols;
		const double hx;
		const double hy;

		std::vector<double> xx;
		std::vector<double> yy;
};


class AnisotropicPoisson2D : public IsotropicPoisson2D {
	public:
		AnisotropicPoisson2D(
			const std::pair<double, double> bottom_left_corner,
			const std::pair<double, double> top_right_corner,
			const int n,
			const std::function<double(double,double)> f,
			const std::function<double(double,double)> boundary,
			const std::function<double(double,double)> c
		);
};


#endif
