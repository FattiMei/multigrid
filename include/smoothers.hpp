#ifndef __SMOOTHERS_HPP__
#define __SMOOTHERS_HPP__


#include "poisson1D.hpp"


namespace Smoother {
class BaseSmoother {
	public:
		BaseSmoother(const Poisson1D &problem);
		virtual double* smooth(const int n, const double b[], double u[]) = 0;


	protected:
		const Update iteration_formula;
};


class Jacobi : public BaseSmoother {
	public:
		Jacobi(const Poisson1D &problem);
		~Jacobi();

		double* smooth(const int n, const double b[], double u[]);


	private:
		double *local;
};


class GSeidel : public BaseSmoother {
	public:
		GSeidel(const Poisson1D &problem) : BaseSmoother(problem) {};
		double* smooth(const int n, const double b[], double u[]);
};


class RedBlack : public BaseSmoother {
	public:
		RedBlack(const Poisson1D &problem) : BaseSmoother(problem) {};
		double* smooth(const int n, const double b[], double u[]);
};
}


#endif
