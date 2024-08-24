#include <iostream>
#include <eigen3/Eigen/Sparse>
#include "poisson.hpp"


int main() {
	IsotropicPoisson1D problem(
		0.0,
		1.0,
		9,
		[](double x){ (void) x; return 0.0; },
		{0.0, 0.0}
	);

	DiscreteOperator *fine   = problem.get_discrete_operator();
	DiscreteOperator *coarse = problem.get_discrete_operator(2);

	Eigen::SparseMatrix<double> A = fine->get_sparse_repr();
	Eigen::SparseMatrix<double> B = coarse->get_sparse_repr();

	std::cout << A << std::endl;
	std::cout << B << std::endl;

	delete fine;
	delete coarse;

	return 0;
}
