#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include "solvers.hpp"
#include "vpointer.hpp"


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
			const Problem*			problem,
			const std::vector<MgOp>		cycle_spec,
			const InitializationStrategy	strategy,
			const UpdateStrategy		smoother,
			RestrictionOperator		restrict = full_weight_restriction,
			ProlongationOperator		prolong  = linear_prolongation
		);
		~MgSolver();

		void step();


	protected:
		void build_level_to_memory_map();

		const int		n;
		const std::vector<MgOp> recipe;
		const UpdateStrategy 	smoother;
		const int		maxlevels;

		const RestrictionOperator  restrict;
		const ProlongationOperator prolong;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> direct_solver;

		std::vector<int>		grid_size;
		std::vector<double*>		grid_solution;
		std::vector<PointerVariant<double>> grid_rhs;
		std::vector<double*>		grid_residual;
		std::vector<DiscreteOperator*>	grid_operator;

		double* solution_memory;
		double* rhs_memory;
		double* residual_memory;
};


int analyze_cycle_recipe(const std::vector<MgOp> &recipe);

// probably needs to be specialized for each problem
std::vector<int> compute_grid_sizes(const int n, const int maxdepth);


namespace MgCycle {
std::vector<MgOp> V(int maxdepth, bool solve = false);
}


#endif
