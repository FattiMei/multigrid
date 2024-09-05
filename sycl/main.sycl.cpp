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
	// SYCL code here...
	cl::sycl::queue q;

	constexpr int n = 257;
	constexpr int depth = 5;
	constexpr int smoothing_steps = 3;
	const auto recipe = MgCycle::V(depth);

	std::vector<int>     grid_dim(depth);
	std::vector<double*> grid_solution(depth);
	std::vector<double*> grid_rhs(depth);
	std::vector<double*> grid_residual(depth);

	// building the memory map
	// initializing the grid memories to 0
	int m = n;
	for (int level = 0; level < depth; ++level) {
		grid_dim[level] = m;

		grid_solution[level] = cl::sycl::malloc_device<double>(m, q);
		grid_rhs     [level] = cl::sycl::malloc_device<double>(m, q);
		grid_residual[level] = cl::sycl::malloc_device<double>(m, q);

		q.fill(grid_solution[level], 0, m * sizeof(double));
		q.fill(grid_rhs     [level], 0, m * sizeof(double));
		q.fill(grid_residual[level], 0, m * sizeof(double));

		m = 1 + (m-1)/2;
	}

	const int coarse_dim = grid_dim.back();

	// host memory for the direct solving
	std::vector<double> coarse_rhs(coarse_dim * coarse_dim);
	std::vector<double> coarse_sol(coarse_dim * coarse_dim);

	const double h = std::pow(2.0, depth) / (n-1.0);
	Eigen::SparseLU<Eigen::SparseMatrix<double>> direct_solver;
	FivePointStencil stencil_op(
		coarse_dim,
		coarse_dim,
		{
			-1.0 / (h*h),
			 4.0 / (h*h),
			-1.0 / (h*h),
			-1.0 / (h*h),
			-1.0 / (h*h)
		}
	);

	direct_solver.compute(stencil_op.get_sparse_repr());

	// building the rhs and transferring to device
	std::vector<double> rhs(n*n);
	for (auto &x : rhs) x = 0.0;

	for (int row = 1; row < n-1; ++row) {
		for (int col = 1; col < n-1; ++col) {
			rhs[row*n+col] = 1.0;
		}
	}
	q.memcpy(grid_rhs[0], rhs.data(), rhs.size() * sizeof(double));
	q.wait();


	// using the fact that the boundary is at zero so no need to iterate on boundary
	
	const auto t1 = std::chrono::high_resolution_clock::now();

	int level = 0;
	for (const auto op : recipe) {
		switch (op) {
			case MgOp::Relax: {
				for (int step = 0; step < smoothing_steps; ++step) {
					cl::sycl::range<2> work_items{static_cast<size_t>(grid_dim[level]-2), static_cast<size_t>(grid_dim[level]-2)};
					const int m = grid_dim[level];

					double* u   = grid_solution[level];
					double* rhs = grid_rhs[level];

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
				}
			} break;

			case MgOp::Restrict: {
				/*
				cl::sycl::range<2> work_items{static_cast<size_t>(grid_dim[level]-2), static_cast<size_t>(grid_dim[level]-2)};
				const int m = grid_dim[level];

				double* u   = grid_solution[level];
				double* rhs = grid_rhs[level];
				double* res = grid_residual[level];
				double* coarse   = grid_residual[level-1];

				q.submit([&](cl::sycl::handler& cgh) {
					cgh.parallel_for<class residual_kernel>(
						work_items,
						[=](cl::sycl::id<2> tid) {
							const int row = tid[0]+1;
							const int col = tid[1]+1;

							const int i = row * m + col;
							const double hx = 1.0 / (m-1.0);

							res[i] = (rhs[i] - u[i-1] + 4.0 * u[i] - u[i+1] - u[i+m] - u[i-m]) / (hx*hx);
						}
					);
				}).wait();

				q.submit([&](cl::sycl::handler& cgh) {
					cgh.parallel_for<class restrict_kernel>(
						work_items,
						[=](cl::sycl::id<2> tid) {
							const int row = tid[0]+1;
							const int col = tid[1]+1;

							const int i = row * m + col;
							const int target_col = 1 + (m-1)/2;

							if (row % 2 == 0 and col % 2 == 0) {
								coarse[(row/2) * target_col + (col/2)] =
									  0.25   * u[i]
									+ 0.125  * u[i - 1]
									+ 0.125  * u[i + 1]
									+ 0.125  * u[i + m]
									+ 0.125  * u[i - m]
									+ 0.0625 * u[i + m + 1]
									+ 0.0625 * u[i + m - 1]
									+ 0.0625 * u[i - m + 1]
									+ 0.0625 * u[i - m - 1];
							}
						}
					);
				}).wait();
				*/

				--level;
			} break;

			case MgOp::Prolong: {
				++level;

			} break;

			case MgOp::DirectSolve: {
				Eigen::Map<Eigen::VectorXd> b  (coarse_rhs.data(), coarse_rhs.size());
				Eigen::Map<Eigen::VectorXd> tmp(coarse_sol.data(), coarse_sol.size());

				q.memcpy(b.data(), grid_rhs[level], b.size());
				q.wait();

				tmp = direct_solver.solve(b);

				q.memcpy(grid_solution[level], coarse_sol.data(), coarse_sol.size());
				q.wait();
			} break;

			case MgOp::IterativeSolve: {

			} break;
		}
	}

	const auto t2 = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double,std::milli> delta = t2 - t1;
	std::cout << delta.count() << "ms" << std::endl;

	return 0;
}
