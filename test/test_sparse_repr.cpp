#include <iostream>
#include <eigen3/Eigen/Sparse>
#include "poisson.hpp"


int main() {
	IsotropicPoisson1D problem(
		0.0,
		1.0,
		5,
		[](double x){ return 0.0; },
		{0.0, 0.0}
	);

	DiscreteOperator *op = problem.get_discrete_operator();
	Eigen::SparseMatrix<double> A = op->get_sparse_repr();

	std::cout << A << std::endl;
	delete op;

	return 0;
}
