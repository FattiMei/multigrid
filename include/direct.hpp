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
