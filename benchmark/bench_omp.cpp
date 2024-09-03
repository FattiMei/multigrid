#include <benchmark/benchmark.h>
#include "multigrid.hpp"
#include "bench_step.hpp"


BENCHMARK(BM_step<5,MgCycle::V,UpdateStrategy::RedBlack,1>)->RangeMultiplier(2)->Range(32, 1024);
BENCHMARK(BM_step<5,MgCycle::V,UpdateStrategy::RedBlack,2>)->RangeMultiplier(2)->Range(32, 1024);
BENCHMARK(BM_step<5,MgCycle::V,UpdateStrategy::RedBlack,4>)->RangeMultiplier(2)->Range(32, 1024);
BENCHMARK(BM_step<5,MgCycle::V,UpdateStrategy::RedBlack,8>)->RangeMultiplier(2)->Range(32, 1024);


BENCHMARK_MAIN();
