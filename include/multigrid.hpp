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
		MgSolver(const Poisson1D &problem, const std::vector<MgOp> recipe, InitializationStrategy strategy);

		void step();


	protected:
		const std::vector<MgOp> recipe;

};


void injective_restriction  (const int n, const double src[], double dest[]);
void full_weight_restriction(const int n, const double src[], double dest[]);
void linear_prolongation    (const int m, const double src[], double dest[]);


bool analyze_cycle_recipe(const std::vector<MgOp> &recipe, int &nlevels);


#endif
