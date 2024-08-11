#ifndef __SMOOTHERS_HPP__
#define __SMOOTHERS_HPP__


#include "poisson1D.hpp"
#include <vector>


namespace Smoother {
class Jacobi {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);


	private:
		// Use std::vector and not a simple pointer to make use of RAII (omitting destructors)
		std::vector<double> local;
};


class GSeidel {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);
};


class RedBlack {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);
};


class BlackRed {
	public:
		void operator () (const Update formula, const int n, const double b[], double u[]);
};
}


#endif
