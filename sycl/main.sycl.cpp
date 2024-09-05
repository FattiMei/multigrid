#include <iostream>
#include <chrono>
#include <CL/sycl.hpp>
#include <eigen3/Eigen/Sparse>
#include "multigrid.hpp"
#include "stencil.hpp"


// @DESIGN: I could include GPU abstraction into the already polymorphic design,
// to minimize code repetitions I need to include allocations and locality as an abstraction:
//   operations on GPU need device memory
//   not all operations can be done on GPU (direct solvers) and so there is the need of transfer
//   stencil need to be generalized to operate on GPU...
//
// All this without really knowing if GPU implementations are beneficial to the overall solve time,
// the sensible approach is to spend as little engineering time on the design side and produce a rough estimation


int main(int argc, char* argv[]) {
	cl::sycl::queue q;

	if (argc != 2) {
		std::cerr << "Usage: ./poisson_sycl <max number of nodes>" << std::endl;
		return 1;
	}

	int maxnodes = std::atoi(argv[1]);

	double *u      = cl::sycl::malloc_device<double>(maxnodes * maxnodes, q);
	double *tmp    = cl::sycl::malloc_device<double>(maxnodes * maxnodes, q);
	double *rhs    = cl::sycl::malloc_device<double>(maxnodes * maxnodes, q);
	double *res    = cl::sycl::malloc_device<double>(maxnodes * maxnodes, q);
	double *coarse = cl::sycl::malloc_device<double>(maxnodes * maxnodes, q);

	std::cout << "n,red-black,jacobi" << std::endl;

	for (int n = 32; n < maxnodes; n *= 2) {
		const int m = n+1;

		std::cout << n << ",";

		{
			const auto t1 = std::chrono::high_resolution_clock::now();


			cl::sycl::range<2> work_items{static_cast<size_t>(m-2), static_cast<size_t>(m-2)};

			// assign only inner nodes to threads
			// use red black ordering but only 50% occupancy
			q.submit([&](cl::sycl::handler& cgh) {
				cgh.parallel_for<class relaxing_kernel>(
					work_items,
					[=](cl::sycl::id<2> tid) {
						const int row = tid[0]+1;
						const int col = tid[1]+1;

						const int i = row * m + col;
						const double hx = 1.0 / (m-1.0);

						if ((row + col) % 2 == 0) {
							u[i] = (hx*hx*rhs[i] - u[i-1] - u[i+1] - u[i+m] - u[i-m]) / 4.0;
						}
					}
				);
			}).wait();

			q.submit([&](cl::sycl::handler& cgh) {
				cgh.parallel_for<class relaxing_kernel>(
					work_items,
					[=](cl::sycl::id<2> tid) {
						const int row = tid[0]+1;
						const int col = tid[1]+1;

						const int i = row * m + col;
						const double hx = 1.0 / (m-1.0);

						if ((row + col) % 2 == 1) {
							u[i] = (hx*hx*rhs[i] - u[i-1] - u[i+1] - u[i+m] - u[i-m]) / 4.0;
						}
					}
				);
			}).wait();

			const auto t2 = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double,std::milli> delta = t2 - t1;
			std::cout << delta.count() << "ms,";
		}
		{
			const auto t1 = std::chrono::high_resolution_clock::now();


			cl::sycl::range<2> work_items{static_cast<size_t>(m-2), static_cast<size_t>(m-2)};

			// assign only inner nodes to threads
			q.submit([&](cl::sycl::handler& cgh) {
				cgh.parallel_for<class relaxing_kernel>(
					work_items,
					[=](cl::sycl::id<2> tid) {
						const int row = tid[0]+1;
						const int col = tid[1]+1;

						const int i = row * m + col;
						const double hx = 1.0 / (m-1.0);

						tmp[i] = (hx*hx*rhs[i] - u[i-1] - u[i+1] - u[i+m] - u[i-m]) / 4.0;
					}
				);
			}).wait();

			const auto t2 = std::chrono::high_resolution_clock::now();
			const std::chrono::duration<double,std::milli> delta = t2 - t1;
			std::cout << delta.count() << "ms" << std::endl;
		}
	}

	return 0;
}
