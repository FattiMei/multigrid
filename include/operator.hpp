#ifndef __OPERATOR_HPP__
#define __OPERATOR_HPP__


#include <eigen3/Eigen/Sparse>
#include <vector>


enum class UpdateStrategy {
	Jacobi,
	GaussSeidel,
	RedBlack,
	SOR
};


class DiscreteOperator {
	public:
		DiscreteOperator(const int nodes) : n(nodes) {};
		virtual ~DiscreteOperator() = default;

		virtual void relax(const double b[], double u[], UpdateStrategy strategy) = 0;
		virtual Eigen::SparseMatrix<double> get_sparse_repr() const = 0;
		virtual void   compute_residual     (const double b[], const double u[], double r[]) const = 0;
		virtual double compute_residual_norm(const double b[], const double u[]) const = 0;


	protected:
		const int n;
};


// @DESIGN: this is not used in the codebase, but provides a possible implementation of more sophisticated, non-stencil based operators
struct SparseOperator {
	public:
		SparseOperator(
			const int n,
			const std::vector<int>    &offset,
			const std::vector<int>    &neighbours,
			const std::vector<double> &weights
		);


	protected:
		const std::vector<int>    &offset;
		const std::vector<int>    &neighbours;
		const std::vector<double> &weights;
};


#endif
