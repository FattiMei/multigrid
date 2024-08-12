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
		~MgSolver();

		void step();


	protected:
		void build_level_to_memory_map();
		const std::vector<MgOp> recipe;
		Smoother::GSeidel smoother;
		int maxlevels;

		std::vector<int>     grid_size;
		std::vector<double*> grid_solution;
		std::vector<double*> grid_rhs;
		std::vector<double*> grid_residual;
		std::vector<Update>  grid_iteration_formula;
		std::vector<Update>  grid_residual_formula;

		double* solution_memory;
		double* rhs_memory;
		double* residual_memory;

};


void injective_restriction  (const int n, const double src[], double dest[]);
void full_weight_restriction(const int n, const double src[], double dest[]);
void linear_prolongation    (const int m, const double src[], double dest[]);


bool analyze_cycle_recipe(const std::vector<MgOp> &recipe, int &nlevels);
bool compute_grid_sizes(const int n, const int maxlevels, std::vector<int> &grid_size);


#endif
