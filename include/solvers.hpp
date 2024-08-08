#ifndef __SOLVERS_HPP__
#define __SOLVERS_HPP__


#include "poisson1D.hpp"


enum class SolverStatus {
	OK,
	MAXIT
};


// a solver that has the ownership of the solution array
class BaseSolver {
	public:
		BaseSolver(const Poisson1D &problem);
		~BaseSolver();


	protected:
		const Poisson1D &problem;
		const int n;
		double *u;
};


// a solver with stepping logic
class IterativeSolver : BaseSolver {
	public:
		IterativeSolver(const Poisson1D &problem);

		virtual void step();
			void solve(const double tol, const int maxiter);


	protected:
		SolverStatus status = SolverStatus::OK;
		int it;
};


#endif
