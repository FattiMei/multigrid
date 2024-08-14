#ifndef __STENCIL_HPP__
#define __STENCIL_HPP__


#include <array>
#include <set>
#include <eigen3/Eigen/Sparse>


struct ThreePointStencilOperator {
	public:
		ThreePointStencilOperator(
			const int problem_size,
			const std::array<double, 3> weights,
			const std::set<int> boundary_set
		) : n(problem_size) , stencil(weights) , boundary(boundary_set) {};


		Eigen::SparseMatrix<double> build_sparse_repr();


	private:
		const int n;
		const std::array<double,3> stencil;
		const std::set<int>        boundary;
};


#endif
