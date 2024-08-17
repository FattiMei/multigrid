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

		DiscreteOperator* get_discrete_operator(const int level = 0) const override;


	private:
		const double h;
};



#endif
