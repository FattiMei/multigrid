#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include "solvers.hpp"


enum class MgOp {
	Relax,
	Restrict,
	Prolong,
	DirectSolve,
	IterativeSolve  // not used in production, only for debugging purpuoses
};


class MgSolver : public IterativeSolver {
	public:
		MgSolver(const Poisson1D &problem);

		void step();


	protected:

};


#endif
