#ifndef __STENCIL_HPP__
#define __STENCIL_HPP__


#include <array>
#include <set>
#include <eigen3/Eigen/Sparse>


// assumes a 1D grid with extremal points at boundary
struct ThreePointStencilOperator {
	public:
		ThreePointStencilOperator(
			const int n,
			const std::array<double, 3> weights
		) : problem_size(n) , stencil(weights) {};


		Eigen::SparseMatrix<double> build_sparse_repr() const;

		const int problem_size;
		const std::array<double,3> stencil;
};


#endif
