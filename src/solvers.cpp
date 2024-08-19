#include "solvers.hpp"
#include "operator.hpp"
#include <iostream>


BaseSolver::BaseSolver(const Problem *p) :
	problem(p),
	op(p->get_discrete_operator()),
	u(p->get_size()),
	rhs(p->get_rhs())
{}


double BaseSolver::get_residual_norm() const {
	return op->compute_residual_norm(rhs, u.data());
}


const std::vector<double>& BaseSolver::get_solution() const {
	return u;
}


BaseSolver::~BaseSolver() {
	delete op;
}


IterativeSolver::IterativeSolver(const Problem *p, InitializationStrategy strategy) : BaseSolver(p) {
	switch (strategy) {
		default: {
			for (auto &x : u) {
				x = 0.0;
			}
		}
	}
}


void IterativeSolver::solve(const double tol, const int maxiter) {
	if (status == SolverStatus::OK) {
		status = SolverStatus::MAXIT;

		for (it = 0; it < maxiter; ++it) {
			this->step();

			if (op->compute_residual_norm(rhs, u.data()) < tol) {
				status = SolverStatus::OK;
				break;
			}
		}
	}
}


SmootherSolver::SmootherSolver(
	const Problem *problem,
	const InitializationStrategy init,
	const UpdateStrategy strategy
) :
	IterativeSolver(problem, init),
	smoother(strategy)
{}


void SmootherSolver::step() {
	op->relax(rhs, u.data(), smoother);
}
