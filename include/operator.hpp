#ifndef __OPERATOR_HPP__
#define __OPERATOR_HPP__


#include <set>
#include <array>
#include <concepts>
#include <eigen3/Eigen/Sparse>


struct SparseOperator {
	public:
		SparseOperator(
			const int n,
			const std::vector<int>    &offset,
			const std::vector<int>    &neighbours,
			const std::vector<double> &weights
		);

		Eigen::SparseMatrix<double> build_sparse_repr();

		const int problem_size;
		const std::vector<int>    &offset;
		const std::vector<int>    &neighbours;
		const std::vector<double> &weights;
};




void compute_residual     (const SparseOperator &A, const double b[], const double u[], double r[]);
void compute_residual_norm(const SparseOperator &A, const double b[], const double u[]);


#endif
