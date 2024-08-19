#ifndef __STENCIL_HPP__
#define __STENCIL_HPP__


#include "operator.hpp"


// assumes a 1D grid with extremal points at boundary
class ThreePointStencil : public DiscreteOperator {
	public:
		ThreePointStencil(const int n, const std::array<double,3> weights);

		void relax(const double b[], double u[], UpdateStrategy strategy) override;
		Eigen::SparseMatrix<double> get_sparse_repr() const override;
		void   compute_residual     (const double b[], const double u[], double r[]) const override;
		double compute_residual_norm(const double b[], const double u[]) const override;


	protected:
		const std::array<double,3> stencil;
		std::vector<double> local;
};


// assumes a 1D grid with extremal points topologically connected
class ThreePointPeriodicStencil : public ThreePointStencil {
	public:
		void relax(const double b[], double u[], UpdateStrategy strategy) override;
		Eigen::SparseMatrix<double> get_sparse_repr() const override;
		void   compute_residual     (const double b[], const double u[], double r[]) const override;
		double compute_residual_norm(const double b[], const double u[]) const override;
};


#endif
