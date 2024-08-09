#ifndef __SMOOTHERS_HPP__
#define __SMOOTHERS_HPP__


#include "poisson1D.hpp"


class Smoother {
	public:
		Smoother(const Poisson1D &problem);
		virtual double* smooth(const int n, const double b[], double u[]) = 0;


	protected:
		const Update iteration_formula;
};


class Jacobi : public Smoother {
	public:
		Jacobi(const Poisson1D &problem);
		~Jacobi();

		double* smooth(const int n, const double b[], double u[]);


	private:
		double *local;
};


class GS : public Smoother {
	public:
		GS(const Poisson1D &problem) : Smoother(problem) {};
		double* smooth(const int n, const double b[], double u[]);
};


class RedBlack : public Smoother {
	public:
		double* smooth(const int n, const double b[], double u[]);
};


#endif
