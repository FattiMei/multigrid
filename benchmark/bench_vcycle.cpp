#include <benchmark/benchmark.h>
#include "symbolic.h"
#include "poisson.hpp"
#include "multigrid.hpp"


template <int DEPTH, UpdateStrategy SMOOTHER>
static void vcycle_step(benchmark::State &state) {
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
		MgCycle::V(DEPTH, 3),
		InitializationStrategy::Zeros,
		SMOOTHER,
		full_weight_restriction_2d,
		linear_prolongation_2d
	);

	for (auto _ : state) {
		solver.step();
	}
}


BENCHMARK(vcycle_step<2, UpdateStrategy::GaussSeidel>)->RangeMultiplier(2)->Range(4, 256);
BENCHMARK(vcycle_step<2, UpdateStrategy::RedBlack>   )->RangeMultiplier(2)->Range(4, 256);
BENCHMARK_MAIN();
