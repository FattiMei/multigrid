#include "operator.hpp"
#include <cmath>


ThreePointStencil::ThreePointStencil(
	const int n,
	const std::array<double,3> weights
) :
	DiscreteOperator(n),
	stencil(weights)
{}


Eigen::SparseMatrix<double> ThreePointStencil::get_sparse_repr() {
	Eigen::SparseMatrix<double> A(n,n);
	A.reserve(Eigen::VectorXi::Constant(n, 3));

	A.coeffRef(0,0) = 1.0;
	for (int i = 1; i < n-1; ++i) {
		A.coeffRef(i,i-1) = stencil[0];
		A.coeffRef(i,i)   = stencil[1];
		A.coeffRef(i,i+1) = stencil[2];
	}
	A.coeffRef(n-1,n-1) = 1.0;

	A.makeCompressed();
	return A;
}


double ThreePointStencil::compute_residual_norm(const double b[], const double u[]) const {
	double acc = 0.0;
	double diff;

	diff = b[0] - u[0];
	acc += diff * diff;

	for (int i = 1; i < n-1; ++i) {
		diff = b[i] - u[i-1] * stencil[0] - u[i] * stencil[1] - u[i+1] * stencil[2];
		acc += diff * diff;
	}

	diff = b[n-1] - u[n-1];
	acc += diff * diff;

	return std::sqrt(acc);
}
