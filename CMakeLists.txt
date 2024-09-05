cmake_minimum_required(VERSION 3.10)
project(Multigrid VERSION 1.0)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# manually insert the path to the adaptiveCpp SYCL compiler (or equivalent SYCL compilers)
set(SYCLCC "/usr/local/bin/acpp")
set(SYCLCCFLAGS "-O3")


# manually apply optimization flags (might experiment with -O3 -ftree-vectorize...)
add_compile_options(-O2)
add_compile_options(-Wall -Wextra -Wpedantic)


find_package(OpenMP REQUIRED)

if (OpenMP_CXX_FOUND)
	add_compile_options(${OpenMP_CXX_FLAGS})
	link_libraries(${OpenMP_CXX_LIBRARIES})
endif()


find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Python3 REQUIRED)

# github.com/google/benchmark for microbenchmarking of C/C++ code and compact result generation
find_package(benchmark REQUIRED)
find_package(AdaptiveCpp CONFIG REQUIRED)


if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
	add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
	add_compile_options(/W4 /WX)
endif()


set(INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include")

set(LIB_SRC_DIR         "${CMAKE_SOURCE_DIR}/src")
set(TEST_SRC_DIR        "${CMAKE_SOURCE_DIR}/test")
set(SYCL_SRC_DIR        "${CMAKE_SOURCE_DIR}/sycl")
set(BENCH_SRC_DIR       "${CMAKE_SOURCE_DIR}/benchmark")
set(EXAMPLES_SRC_DIR    "${CMAKE_SOURCE_DIR}/examples")
set(ARTIFACTS_SRC_DIR   "${CMAKE_SOURCE_DIR}/artifacts")

file(GLOB LIB_SRC       "${LIB_SRC_DIR}/*.cpp")
file(GLOB TEST_SRC      "${TEST_SRC_DIR}/test_*.cpp")
file(GLOB SYCL_SRC      "${SYCL_SRC_DIR}/*.sycl.cpp")
file(GLOB BENCH_SRC     "${BENCH_SRC_DIR}/*.cpp")
file(GLOB EXAMPLES_SRC  "${EXAMPLES_SRC_DIR}/example_*.cpp")
file(GLOB ARTIFACTS_SRC "${ARTIFACTS_SRC_DIR}/*")


add_library(MgLib STATIC ${LIB_SRC})
include_directories(MgLib PRIVATE ${INCLUDE_DIRS})


add_library(Exact STATIC ${ARTIFACTS_SRC})
include_directories(Exact PRIVATE ${ARTIFACTS_SRC_DIR})


enable_testing()


foreach(TEST_SOURCE ${TEST_SRC})
	get_filename_component(EXE_NAME ${TEST_SOURCE} NAME_WE)

	add_executable(${EXE_NAME} ${TEST_SOURCE})
	target_link_libraries(${EXE_NAME} PRIVATE MgLib)

	add_test(NAME ${EXE_NAME} COMMAND ${EXE_NAME})
endforeach()


foreach(BENCH_SOURCE ${BENCH_SRC})
	get_filename_component(EXE_NAME ${BENCH_SOURCE} NAME_WE)

	add_executable(${EXE_NAME} ${BENCH_SOURCE})
	target_link_libraries(${EXE_NAME} PRIVATE MgLib Exact benchmark::benchmark)
endforeach()


foreach(EXAMPLE_SOURCE ${EXAMPLES_SRC})
	get_filename_component(EXE_NAME ${EXAMPLE_SOURCE} NAME_WE)

	add_executable(${EXE_NAME} ${EXAMPLE_SOURCE})
	target_include_directories(${EXE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/artifacts)
	target_link_libraries(${EXE_NAME} PRIVATE MgLib)
	target_link_libraries(${EXE_NAME} PRIVATE Exact)
endforeach()


# produce the exact solution and forcing term to be used in convergence tests
add_custom_command(
	OUTPUT ${ARTIFACTS_SRC}
	DEPENDS ${CMAKE_SOURCE_DIR}/analysis/poisson_exact_solution.py
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/analysis/poisson_exact_solution.py ${CMAKE_SOURCE_DIR}/artifacts/symbolic.cpp ${CMAKE_SOURCE_DIR}/artifacts/symbolic.h
	COMMENT "Generating C functions from symbolic expressions"
)

add_custom_target(
	run_symbolic_generation ALL
	DEPENDS ${CMAKE_SOURCE_DIR}/artifacts/symbolic.cpp
	DEPENDS ${CMAKE_SOURCE_DIR}/artifacts/symbolic.h
)


# custom rules for calling experiments, saving the results and calling the proper python script to plot them
add_custom_target(
	circuit
	DEPENDS example_forced_circuit
	COMMAND example_forced_circuit > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/solution_plot_1D.py out
)


add_custom_target(
	multilevel
	DEPENDS example_multilevel_performance
	COMMAND example_multilevel_performance > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/convergence_history.py out
)


add_custom_target(
	cycle
	DEPENDS example_cycle_performance
	COMMAND example_cycle_performance > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/convergence_history.py out
)


add_custom_target(
	convergence_1d
	DEPENDS example_convergence_1d
	COMMAND example_convergence_1d > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/convergence_order.py out
)


add_custom_target(
	convergence_2d
	DEPENDS example_convergence_2d
	COMMAND example_convergence_2d > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/convergence_order.py out
)


add_custom_target(
	smoother
	DEPENDS example_smoother_comparison
	COMMAND example_smoother_comparison > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/convergence_history.py out
)


# custom rules for benchmarking
add_custom_target(
	benchmark_omp
	DEPENDS bench_omp
	COMMAND ./bench_omp --benchmark_out=omp.json --benchmark_out_format=json
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/speedup.py omp.json
)


add_custom_target(
	benchmark_vcycle
	DEPENDS bench_vcycle
	COMMAND ./bench_vcycle --benchmark_out=vcycle.json --benchmark_out_format=json
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/speedup.py vcycle.json
)


add_custom_target(
	benchmark_cycle
	DEPENDS bench_cycle
	COMMAND ./bench_cycle --benchmark_out=cycle.json --benchmark_out_format=json
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/speedup.py cycle.json
)


add_custom_target(
	benchmark_fused
	DEPENDS bench_fused
	COMMAND ./bench_fused --benchmark_out=fused.json --benchmark_out_format=json
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/speedup.py fused.json
)


add_custom_target(
	poisson_sycl ALL
	COMMAND ${SYCLCC} ${SYCLCCFLAGS} -I ${CMAKE_SOURCE_DIR}/include -I ${CMAKE_SOURCE_DIR}/sycl -o poisson_sycl ${SYCL_SRC}
)
