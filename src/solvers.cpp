#include "solvers.hpp"
#include "operator.hpp"
#include <iostream>


/*
BaseSolver::BaseSolver(const Poisson1D &problem_) : problem(problem_), n(problem.get_problem_size()), rhs(problem.get_rhs()) {
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
*/


BaseSolver::BaseSolver(const Problem *p) :
	problem(p),
	op(p->get_discrete_operator()),
	u(p->get_size()),
	rhs(p->get_rhs())
{}


double BaseSolver::get_residual_norm() const {
	return op->compute_residual_norm(rhs, u.data());
}


BaseSolver::~BaseSolver() {
	delete op;
}


EigenDirectSolver::EigenDirectSolver(const Problem *p) :
	problem(p),
	op(problem->get_discrete_operator()),
	u(problem->get_size()),
	rhs(problem->get_rhs())
{
	A = op->get_sparse_repr();
	solver.compute(A);

	if (solver.info() != Eigen::Success) {
		std::cerr << "EigenDirectSolver has failed to factorize the operator" << std::endl;
	}
}


EigenDirectSolver::~EigenDirectSolver() {
	delete op;
}


void EigenDirectSolver::solve() {
	Eigen::Map<const Eigen::VectorXd> b(rhs, problem->get_size());
	Eigen::Map<Eigen::VectorXd> tmp(u.data(), u.size());

	tmp = solver.solve(b);
}


double EigenDirectSolver::get_residual_norm() const {
	return op->compute_residual_norm(rhs, u.data());
}
