#ifndef __DIRECT_HPP__
#define __DIRECT_HPP__


#include "solvers.hpp"
#include <iostream>
#include <eigen3/Eigen/Sparse>


template <class EigenSparseSolver>
class DirectSolver : public BaseSolver {
	public:
		DirectSolver(const Problem *problem) : BaseSolver(problem) {
			solver.compute(op->get_sparse_repr());
			if (solver.info() != Eigen::Success) {
				std::cerr << "EigenDirectSolver has failed to factorize the operator" << std::endl;
			}
		}

		void solve() {
			Eigen::Map<const Eigen::VectorXd> b(rhs, problem->get_size());
			Eigen::Map<Eigen::VectorXd> tmp(u.data(), u.size());

			tmp = solver.solve(b);
		}


	private:
		EigenSparseSolver solver;
};


template <class EigenSparseSolver>
class SymmetricSolver : public BaseSolver {
	public:
		SymmetricSolver(const Problem *problem) : BaseSolver(problem), n(problem->get_size()), h(problem->get_step()) {
			const Eigen::SparseMatrix<double> A = op->get_sparse_repr();
			const Eigen::SparseMatrix<double> inner = A.block(1,1,A.rows()-2,A.cols()-2);

			solver.compute(inner);
			if (solver.info() != Eigen::Success) {
				std::cerr << "EigenDirectSolver has failed to factorize the operator" << std::endl;
			}
		}

		void solve() {
			Eigen::Map<const Eigen::VectorXd> b(rhs+1, problem->get_size()-2);
			Eigen::Map<Eigen::VectorXd> tmp(u.data()+1, u.size()-2);

			tmp = solver.solve(b);
		}


	private:
		const double n;
		const double h;
		EigenSparseSolver solver;
};


#endif
