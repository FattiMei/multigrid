#include <benchmark/benchmark.h>
#include "multigrid.hpp"
#include "bench_step.hpp"


BENCHMARK(BM_step      <5,MgCycle::V,UpdateStrategy::RedBlack,4>)->RangeMultiplier(2)->Range(32,2048);
BENCHMARK(BM_step_fused<5,MgCycle::VF,UpdateStrategy::RedBlack,4>)->RangeMultiplier(2)->Range(32,2048);
BENCHMARK_MAIN();
