#ifndef __SMOOTHERS_HPP__
#define __SMOOTHERS_HPP__


#include "poisson1D.hpp"


namespace Smoother {
class BaseSmoother {
	public:
		BaseSmoother(const int n, const Update iteration_formula);
		virtual void smooth(const double b[], double u[]) = 0;


	protected:
		const Update iteration_formula;
		const int n;
};


class Jacobi : public BaseSmoother {
	public:
		Jacobi(const int n, const Update iteration_formula);
		~Jacobi();

		void smooth(const double b[], double u[]);


	private:
		double *local;
};


class GSeidel : public BaseSmoother {
	public:
		GSeidel(const int n, const Update iteration_formula) : BaseSmoother(n, iteration_formula) {};
		void smooth(const double b[], double u[]);
};


class RedBlack : public BaseSmoother {
	public:
		RedBlack(const int n, const Update iteration_formula) : BaseSmoother(n, iteration_formula) {};
		void smooth(const double b[], double u[]);
};
}


#endif
