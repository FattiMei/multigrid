#ifndef __MULTIGRID_HPP__
#define __MULTIGRID_HPP__


#include "solvers.hpp"
#include "vpointer.hpp"
#include <iostream>


enum class MgOp {
	Relax,
	Restrict,
	Prolong,
	DirectSolve,
	IterativeSolve  // not used in production, only for debugging purpuoses
};


std::ostream& operator<<(std::ostream& os, const MgOp op);


using RestrictionOperator  = std::function<void(const std::pair<int,int>, const double*, double*)>;
using ProlongationOperator = std::function<void(const std::pair<int,int>, const double*, double*)>;


void injection_restriction_1d  (const std::pair<int,int> dim, const double src[], double dest[]);
void full_weight_restriction_1d(const std::pair<int,int> dim, const double src[], double dest[]);
void linear_prolongation_1d    (const std::pair<int,int> dim, const double src[], double dest[]);

void injection_restriction_2d  (const std::pair<int,int> dim, const double src[], double dest[]);
void full_weight_restriction_2d(const std::pair<int,int> dim, const double src[], double dest[]);
void linear_prolongation_2d    (const std::pair<int,int> dim, const double src[], double dest[]);

void linear_prolongation_and_correction_2d(const std::pair<int,int> dim, const double src[], double dest[]);


class MgSolver : public IterativeSolver {
	public:
		MgSolver(
			const Problem*			problem,
			const std::vector<MgOp>		cycle_spec,
			const InitializationStrategy	strategy,
			const UpdateStrategy		smoother,
			RestrictionOperator		restrict,
			ProlongationOperator		prolong
		);
		~MgSolver();

		void step() override;


	protected:
		void build_level_to_memory_map();

		const int		n;
		const std::vector<MgOp> recipe;
		const UpdateStrategy 	smoother;
		const int		maxlevels;

		const RestrictionOperator  restrict;
		const ProlongationOperator prolong;

		Eigen::SparseLU<Eigen::SparseMatrix<double>> direct_solver;

		std::vector<std::pair<int,int>>		grid_dim;
		std::vector<int>			grid_size;
		std::vector<double*>			grid_solution;
		std::vector<PointerVariant<double>>	grid_rhs;
		std::vector<double*>			grid_residual;
		std::vector<DiscreteOperator*>		grid_operator;

		double* solution_memory;
		double* rhs_memory;
		double* residual_memory;
};


class MgFusedSolver : public MgSolver {
	public:
		MgFusedSolver(
			const Problem*			problem,
			const std::vector<MgOp>		cycle_spec,
			const InitializationStrategy	strategy,
			const UpdateStrategy		smoother,
			RestrictionOperator		restrict,
			ProlongationOperator		prolong_and_correct
		) : MgSolver(problem, cycle_spec, strategy, smoother, restrict, prolong_and_correct) {};

		void step() override;
};


int analyze_cycle_recipe(const std::vector<MgOp> &recipe);

std::vector<std::pair<int,int>> compute_grid_dim(const Problem* problem, const int maxdepth);
std::vector<int> compute_grid_size(const std::vector<std::pair<int,int>>& sizes);


namespace MgCycle {
	// @DESIGN: levels is one indexed so that levels=2 for a multigrid with two grids (main and coarse). Sacrifice ease of implementation for usability
	std::vector<MgOp> V(const int levels, const int smoothing_steps = 1, const bool solve = true);
	std::vector<MgOp> F(const int levels, const int smoothing_steps = 1, const bool solve = true);
	std::vector<MgOp> W(const int levels, const int smoothing_steps = 1, const bool solve = true);
}


#endif
