#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include "solvers.hpp"
#include "pointer_variant.hpp"


enum class MgOp {
	Relax,
	Restrict,
	Prolong,
	DirectSolve,
	IterativeSolve  // not used in production, only for debugging purpuoses
};


using RestrictionOperator  = std::function<void(const int, const double*, double*)>;
using ProlongationOperator = std::function<void(const int, const double*, double*)>;


void injective_restriction  (const int n, const double src[], double dest[]);
void full_weight_restriction(const int n, const double src[], double dest[]);
void linear_prolongation    (const int m, const double src[], double dest[]);


class MgSolver : public IterativeSolver {
	public:
		MgSolver(
			const Poisson1D &problem,
			const std::vector<MgOp> recipe,
			InitializationStrategy strategy,
			RestrictionOperator restrict = full_weight_restriction,
			ProlongationOperator prolong = linear_prolongation
		);
		~MgSolver();

		void step();


	protected:
		void build_level_to_memory_map();
		const std::vector<MgOp> recipe;
		Smoother::GSeidel smoother;
		int maxlevels;

		RestrictionOperator  restrict;
		ProlongationOperator prolong;

		std::vector<int>     grid_size;
		std::vector<double*> grid_solution;
		std::vector<PointerVariant<double>> grid_rhs;
		std::vector<double*> grid_residual;
		std::vector<Update>  grid_iteration_formula;
		std::vector<Update>  grid_residual_formula;

		double* solution_memory;
		double* rhs_memory;
		double* residual_memory;
};


bool analyze_cycle_recipe(const std::vector<MgOp> &recipe, int &nlevels);
bool compute_grid_sizes(const int n, const int maxlevels, std::vector<int> &grid_size);


namespace MgCycle {
std::vector<MgOp> V(int maxdepth);
}



#endif
