#ifndef __SOLVERS_HPP__
#define __SOLVERS_HPP__


#include "problem.hpp"


enum class SolverStatus {
	OK,
	MAXIT
};


class BaseSolver {
	public:
		BaseSolver(const Problem *problem);
		virtual ~BaseSolver();

		double get_residual_norm() const;


	protected:
		const Problem* problem;
		DiscreteOperator* op;

		std::vector<double> u;
		const double* rhs;
};


// a solver with stepping logic
class IterativeSolver : public BaseSolver {
	public:
		IterativeSolver(
			const Problem *problem,
			const InitializationStrategy strategy
		);
		virtual	void step() = 0;
		void solve(const double tol, const int maxiter);


	protected:
		SolverStatus status = SolverStatus::OK;
		int it;
};


class SmootherSolver : public IterativeSolver {
	public:
		SmootherSolver(
			const Problem *problem,
			const InitializationStrategy init,
			const UpdateStrategy strategy
		);
		void step() override;


	private:
		const UpdateStrategy smoother;
};


#endif
