#ifndef __SMOOTHERS_HPP__
#define __SMOOTHERS_HPP__


#include "poisson1D.hpp"
#include "operator.hpp"
#include "stencil.hpp"
#include <vector>


namespace Smoother {
class Jacobi {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);
		void operator () (const SparseOperator &A           , const double b[], double u[]);
		void operator () (const ThreePointStencilOperator &A, const double b[], double u[]);


	private:
		std::vector<double> local;
};


class GSeidel {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);
		void operator () (const SparseOperator &A           , const double b[], double u[]);
		void operator () (const ThreePointStencilOperator &A, const double b[], double u[]);
};


class RedBlack {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);
		void operator () (const SparseOperator &A           , const double b[], double u[]);
		void operator () (const ThreePointStencilOperator &A, const double b[], double u[]);
};


class BlackRed {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);
		void operator () (const SparseOperator &A           , const double b[], double u[]);
		void operator () (const ThreePointStencilOperator &A, const double b[], double u[]);
};
}


#endif
