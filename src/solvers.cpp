#include "solvers.hpp"


BaseSolver::BaseSolver(const Poisson1D &problem_) : problem(problem_), n(problem.get_problem_size()) {
	u = new double[n];
}


BaseSolver::~BaseSolver() {
	delete[] u;
}


IterativeSolver::IterativeSolver(const Poisson1D &problem, InitializationStrategy strategy) : BaseSolver(problem) {
	problem.set_initial_approximation(u, strategy);
}


double IterativeSolver::get_residual_norm() {
	return problem.get_residual_norm(u);
}


void IterativeSolver::solve(const double tol, const int maxiter) {
	if (status == SolverStatus::OK) {
		status = SolverStatus::MAXIT;

		for (it = 0; it < maxiter; ++it) {
			this->step();

			if (problem.get_residual_norm(u) < tol) {
				status = SolverStatus::OK;
				break;
			}
		}
	}
}


void GsSolver::step() {
	smoother.smooth(n, problem.get_rhs(), u);
}
