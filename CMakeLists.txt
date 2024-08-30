cmake_minimum_required(VERSION 3.10)
project(Multigrid VERSION 1.0)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# do a manual optimization
add_compile_options(-O2)
add_compile_options(-Wall -Wextra -Wpedantic)


find_package(OpenMP REQUIRED)

if (OpenMP_CXX_FOUND)
	add_compile_options(${OpenMP_CXX_FLAGS})
	link_libraries(${OpenMP_CXX_LIBRARIES})
endif()


find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if (Eigen3_CXX_FOUND)
	add_compile_options(${Eigen3_CXX_FLAGS})
endif()


find_package(Python3 REQUIRED)


if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
	add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
	add_compile_options(/W4 /WX)
endif()


set(INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include")

set(LIB_SRC_DIR      "${CMAKE_SOURCE_DIR}/src")
set(TEST_SRC_DIR     "${CMAKE_SOURCE_DIR}/test")
set(EXAMPLES_SRC_DIR "${CMAKE_SOURCE_DIR}/examples")

file(GLOB LIB_SRC      "${LIB_SRC_DIR}/*.cpp")
file(GLOB TEST_SRC     "${TEST_SRC_DIR}/test_*.cpp")
file(GLOB EXAMPLES_SRC "${EXAMPLES_SRC_DIR}/example_*.cpp")


add_library(MgLib STATIC ${LIB_SRC})
target_include_directories(MgLib PRIVATE ${INCLUDE_DIRS})


enable_testing()


foreach(TEST_SOURCE ${TEST_SRC})
	get_filename_component(EXE_NAME ${TEST_SOURCE} NAME_WE)

	add_executable(${EXE_NAME} ${TEST_SOURCE})
	target_link_libraries(${EXE_NAME} PRIVATE MgLib)
	target_include_directories(${EXE_NAME} PRIVATE ${INCLUDE_DIRS})

	add_test(NAME ${EXE_NAME} COMMAND ${EXE_NAME})
endforeach()


foreach(EXAMPLE_SOURCE ${EXAMPLES_SRC})
	get_filename_component(EXE_NAME ${EXAMPLE_SOURCE} NAME_WE)

	add_executable(${EXE_NAME} ${EXAMPLE_SOURCE})
	target_link_libraries(${EXE_NAME} PRIVATE MgLib)
	target_include_directories(${EXE_NAME} PRIVATE ${INCLUDE_DIRS})
endforeach()


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
	convergence
	DEPENDS example_convergence
	COMMAND example_convergence > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/convergence_order.py out
)


add_custom_target(
	smoother
	DEPENDS example_smoother_comparison
	COMMAND example_smoother_comparison > out
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/plot/convergence_history.py out
)
