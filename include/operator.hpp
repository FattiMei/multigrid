#ifndef __OPERATOR_HPP__
#define __OPERATOR_HPP__


#include <eigen3/Eigen/Sparse>


class DiscreteOperator {
	public:
		DiscreteOperator(const int nodes) : n(nodes) {};
		virtual ~DiscreteOperator() = default;

		virtual Eigen::SparseMatrix<double> get_sparse_repr() const = 0;
		virtual double compute_residual_norm(const double b[], const double u[]) const = 0;


	protected:
		const int n;
};


// assumes a 1D grid with extremal points at boundary
class ThreePointStencil : public DiscreteOperator {
	public:
		ThreePointStencil(const int n, const std::array<double,3> weights);

		Eigen::SparseMatrix<double> get_sparse_repr() const override;
		double compute_residual_norm(const double b[], const double u[]) const override;


	private:
		const std::array<double,3> stencil;
};


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
