#ifndef __BENCH_STEP_HPP
#define __BENCH_STEP_HPP


#include "symbolic.h"
#include "poisson.hpp"


template <int DEPTH, auto CYCLE, auto SMOOTHER, int OMP_NUM_THREADS=1>
static void BM_step(benchmark::State &state) {
	omp_set_num_threads(OMP_NUM_THREADS);
	const int n = state.range(0) + 1;

	IsotropicPoisson2D problem(
		{0.0, 0.0},
		{1.0, 1.0},
		n,
		forcing_term_2d,
		solution_2d
	);

	MgSolver solver(
		&problem,
		CYCLE(DEPTH, 3, true),
		InitializationStrategy::Zeros,
		SMOOTHER,
		full_weight_restriction_2d,
		linear_prolongation_2d
	);

	for (auto _ : state) {
		solver.step();
	}
}


#endif
