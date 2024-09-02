#include <benchmark/benchmark.h>
#include "multigrid.hpp"
#include "bench_step.hpp"


BENCHMARK(BM_step<2,MgCycle::V,UpdateStrategy::GaussSeidel>)->RangeMultiplier(2)->Range(4, 256);
BENCHMARK(BM_step<2,MgCycle::V,UpdateStrategy::RedBlack>   )->RangeMultiplier(2)->Range(4, 256);
BENCHMARK(BM_step<2,MgCycle::V,UpdateStrategy::SOR>        )->RangeMultiplier(2)->Range(4, 256);
BENCHMARK_MAIN();
