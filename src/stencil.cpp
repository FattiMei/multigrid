#include "stencil.hpp"


Eigen::SparseMatrix<double> ThreePointStencilOperator::build_sparse_repr() const {
	const int n = problem_size;
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
