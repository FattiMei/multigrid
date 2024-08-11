#ifndef __SOLVERS_HPP__
#define __SOLVERS_HPP__


#include "poisson1D.hpp"
#include "smoothers.hpp"


enum class SolverStatus {
	OK,
	MAXIT
};


// a solver that has the ownership of the solution array
class BaseSolver {
	public:
		BaseSolver(const Poisson1D &problem);
		virtual ~BaseSolver();

		void	set_rhs(const double *rhs);


	protected:
		const Poisson1D &problem;
		const int n;
		double *u;
		const double *rhs;
};


// a solver with stepping logic
class IterativeSolver : public BaseSolver {
	public:
		IterativeSolver(const Poisson1D &problem, InitializationStrategy strategy);

			double	get_residual_norm();
		virtual	void	step() = 0;
			void	solve(const double tol, const int maxiter);


	protected:
		SolverStatus status = SolverStatus::OK;
		int it;
};


template <class Smoother>
class SmootherSolver : public IterativeSolver {
	public:
		SmootherSolver(const Poisson1D &problem, InitializationStrategy strategy) : IterativeSolver(problem, strategy), formula(problem.get_iteration_formula()) {};
		void step() {
			smoother(formula, n, rhs, u);
		};


	private:
		const Update formula;
		Smoother smoother;
};


#endif
