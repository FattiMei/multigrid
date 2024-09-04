#include <benchmark/benchmark.h>
#include "multigrid.hpp"
#include "bench_step.hpp"


BENCHMARK(BM_step<5,MgCycle::V,UpdateStrategy::GaussSeidel>)->RangeMultiplier(2)->Range(32, 1024);
BENCHMARK(BM_step<5,MgCycle::F,UpdateStrategy::GaussSeidel>)->RangeMultiplier(2)->Range(32, 1024);
BENCHMARK(BM_step<5,MgCycle::W,UpdateStrategy::GaussSeidel>)->RangeMultiplier(2)->Range(32, 1024);

BENCHMARK_MAIN();
