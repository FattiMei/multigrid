#include "stencil.hpp"


Eigen::SparseMatrix<double> StencilOperator::build_sparse_repr() {
	Eigen::SparseMatrix<double> A(n,n);
	A.reserve(Eigen::VectorXi::Constant(n, 3));

	for (int i = 0; i < n; ++i) {
		if (boundary.contains(i)) {
			A.coeffRef(i,i) = 1.0;
		}
		else {
			A.coeffRef(i, (i-1) % n) = stencil[0];
			A.coeffRef(i, i)         = stencil[1];
			A.coeffRef(i, (i+1) % n) = stencil[2];
		}
	}

	A.makeCompressed();
	return A;
}
