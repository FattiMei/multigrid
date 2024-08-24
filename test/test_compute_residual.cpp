#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "poisson.hpp"
#include "utils.hpp"


double eigen_compute_residual_norm(const Eigen::SparseMatrix<double> &A, const double rhs[], const double u[]) {
	const int n = A.rows();
	const Eigen::Map<const Eigen::VectorXd> b(rhs, n);
	const Eigen::Map<const Eigen::VectorXd> x(u, n);

	return (b - A * x).norm();
}


int main() {
	constexpr int n = 1000;
	std::vector<double> solution = generate_random_vector(n, 0.0, 1.0);

	IsotropicPoisson1D problem(
		0.0,
		1.0,
		n,
		[](double x){ (void) x; return 0.0; },
		{0.0, 0.0}
	);

	DiscreteOperator *op = problem.get_discrete_operator();
	Eigen::SparseMatrix<double> A = op->get_sparse_repr();

	const double reference   = eigen_compute_residual_norm(A, problem.get_rhs(), solution.data());
	const double alternative = op->compute_residual_norm(problem.get_rhs(), solution.data());

	assert(std::abs((reference - alternative) / reference) < 1e-10);

	delete op;
	return 0;
}
