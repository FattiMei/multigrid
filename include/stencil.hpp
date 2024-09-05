#ifndef __STENCIL_HPP__
#define __STENCIL_HPP__


#include "operator.hpp"


// assumes a 1D grid with extremal points at boundary
class NaiveThreePointStencil : public DiscreteOperator {
	public:
		NaiveThreePointStencil(const int n, const std::array<double,3> weights);

		void relax(const double b[], double u[], UpdateStrategy strategy) override;
		Eigen::SparseMatrix<double> get_sparse_repr() const override;
		void   compute_residual     (const double b[], const double u[], double r[]) const override;
		double compute_residual_norm(const double b[], const double u[]) const override;


	protected:
		const std::array<double,3> stencil;
		std::vector<double> local;
};


// assumes a 1D grid with extremal points at boundary, employs a formula to minimize round-off
class ThreePointStencil : public DiscreteOperator {
	public:
		ThreePointStencil(const int n, const double h, const std::array<double,3> weights);

		void relax(const double b[], double u[], UpdateStrategy strategy) override;
		Eigen::SparseMatrix<double> get_sparse_repr() const override;
		void   compute_residual     (const double b[], const double u[], double r[]) const override;
		double compute_residual_norm(const double b[], const double u[]) const override;


	protected:
		const std::array<double,3> stencil;
		const double h;
		std::vector<double> local;
};


// assumes a 1D grid with extremal points topologically connected with ghost node
class ThreePointPeriodicStencil : public NaiveThreePointStencil {
	public:
		ThreePointPeriodicStencil(const int n, const std::array<double,3> weights) : NaiveThreePointStencil(n, weights) {};

		void relax(const double b[], double u[], UpdateStrategy strategy) override;
		Eigen::SparseMatrix<double> get_sparse_repr() const override;
		void   compute_residual     (const double b[], const double u[], double r[]) const override;
		double compute_residual_norm(const double b[], const double u[]) const override;
};


// assumes a 2D square grid with extremal points at boundary
class FivePointStencil : public DiscreteOperator {
	public:
		FivePointStencil(const int rows, const int cols, const std::array<double,5> weights);

		void relax(const double b[], double u[], UpdateStrategy strategy) override;
		Eigen::SparseMatrix<double> get_sparse_repr() const override;
		void   compute_residual     (const double b[], const double u[], double r[]) const override;
		double compute_residual_norm(const double b[], const double u[]) const override;

		void compute_residual_and_restrict(const std::pair<int,int> dim, const double b[], const double u[], double dest[]);


	protected:
		const std::array<double,5> stencil;
		const int rows;
		const int cols;
		std::vector<double> local;
};


#endif
