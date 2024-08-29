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
