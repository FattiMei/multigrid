#include <benchmark/benchmark.h>
#include "symbolic.h"
#include "poisson.hpp"
#include "multigrid.hpp"


template <UpdateStrategy SMOOTHER>
static void smoother(benchmark::State &state) {
	omp_set_num_threads(4);
	const int n = state.range(0) + 1;

	IsotropicPoisson2D problem(
		{0.0, 0.0},
		{1.0, 1.0},
		n,
		forcing_term_2d,
		solution_2d
	);

	SmootherSolver solver(
		&problem,
		InitializationStrategy::Zeros,
		SMOOTHER
	);

	for (auto _ : state) {
		solver.step();
	}
}


BENCHMARK(smoother<UpdateStrategy::GaussSeidel>)->RangeMultiplier(2)->Range(32,8192);
BENCHMARK(smoother<UpdateStrategy::RedBlack>   )->RangeMultiplier(2)->Range(32,8192);
BENCHMARK(smoother<UpdateStrategy::SOR>        )->RangeMultiplier(2)->Range(32,8192);

BENCHMARK_MAIN();
