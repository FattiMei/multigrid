# multigrid
This project has started as a multigrid implementation of the Poisson problem for regular 2D grids. I'm extending it to include convolutions as the basic building block of the solver algorithm

## Mathematical description
The Poisson problem is a partial differential equation of the form:
$$ \nabla \cdot (\nabla u) = f $$

we consider the isotropic case on a rectangular domain:
$$ \nabla^2 u = f $$

the resulting problem has very rigid assumptions but it could still be applied for solving the pressure predictor step in [incompressible flows simulations](https://github.com/FattiMei/mpi-incompressible-fluid).

### Boundary conditions
Here is a list of boundary conditions that I plan to support
  * homogeneous Dirichlet (aka zero at the border)
  * non-homogeneous Dirichlet

The latter can always be transformed in the homogeneous case using the following transformation:
  1. let $g(x,y)$ be the boundary function
  2. produce a $w(x,y)$ that interpolate $g$ on the problem domain ([bilinear interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation) should be good)
  3. solve the Poisson problem for $\tilde{u} - w$ which now is zero at the border

The final transformation yields a new PDE:
$$ \nabla^2 \tilde{u} = f - \nabla^2 w $$

## Solver zoo
After a finite difference discretization we are left with a linear system $Ax=b$ to solve:
  * **direct solver** uses a form of gaussian elimination. Very good in the 1D case, can't scale well in the 2D and 3D cases due to fill-in phenomenon
  * **fast Poisson solver** diagonalizes the matrix $A$, exploiting its structure
  * black box **iterative solvers** like gradient and conjugate gradient
  * **multigrid** combines an iterative part and a direct solver

In theory the fast Poisson solver has the best algorithmic complexity, but some exploration must be done to measure the constants. It diagonalizes the matrix $A$
$$ A = S^{-1} \Lambda S $$

the columns of $S$ are the eigenvectors of $A$ i.e the eigenstates of the dicrete laplacian operator. With this factorization the system solution becomes:
$$ x = S^{-1} \Lambda^{-1} S b $$

Given the structure of $A$ it's eigenvectors (when homogeneous Dirichlet boundary conditions are enforced) are the base of sines. The *Discrete Sine Transform* and it's inverse are an efficient way of computing the *action* of $S$ and $S^{-1}$, thus giving a speedup over traditional methods. **REMARK:** such method relies on assumptions about the regularity of the problem both in the operator and in the mesh.

## Possible edge for multigrid methods
The multigrid method works by solving the error equation in coarser grids. The coarsest problem has to be solved with high accuracy to ensure the proper convergence behaviour so direct methods must be used. 

There is a balancing between performing FFTs on large grids (thus solving exactly the problem) and using multigrid to reduce the dimension of the grid. I estimate that for large enough inputs the multigrid approach will be the optimal one. This will hold if the smoothing-restriction-prolongation operators have a low impact on the iteration time.

## Describing multigrid primitives as convolutions
When discretizing an operator $\mathcal{L}$ using the finite difference method on a regular grid, the resulting matrix $A$ acts on the solution vector $u$ as a convolution.

Convolution at the boundary nodes has to be treated as a special case. For Dirichlet conditions the boundary nodes are fixed, so they can be omitted in the solution tensor but their neighbours still need the information about the boundary value. One would have to either:
  * implement special cases for the convolution at the border
  * work with the complete tensor and iterate only on the internal nodes

those solutions can't be easily integrated with the convolutions libraries. Fortunately as we saw above every non-homogeneous Dirichlet problem can be transformed into an homogeneous one. This means that we can work with the tensor of internal points and use *zero padding* in the convolution functions.

### Benefits of convolutions
Recent developments in deep learning, GPUs and accelerators have improved the performance of convolutions. Some chips even have a NPU (neural processing units) that implements natively convolutions (although with at most fp16 precision). This means that there are plenty of fast `conv2d` implementations for every type of hardware, we can exploit the engineering effort for solving scientific problems, regardless of the target hardware (vectorized CPU, GPU, NPU, FPGA...)

Here is a list of candidate libraries:
  * pyTorch
  * TensorFlow
  * jax
  * CUDA related


## Build
Has the following dependencies:
  * *C++20 compiler*
  * *OpenMP:* for cpu parallelization of solvers
  * *Eigen*: for direct solving of sparse linear systems (>=3.3)
  * [google benchmark:](https://github.com/google/benchmark) for measuring solvers performance (last version should work fine)
  * *python 3 + numpy + matplotlib + pandas*: for plotting experimental results
  * [sympy](https://www.sympy.org): for producing symbolically the exact solution for convergence tests

Compile the examples and tests as a standard cmake project:
```(bash)
mkdir build
cd build
cmake ..
make
```

The `SYCL` branch has a SYCL compiler dependency, further compilation instructions are in `CMakeLists.txt`.


## Overview
The project implements a geometric multigrid solver for 1D and 2D square grids. It uses C++ runtime polimorphism to capture the variability of PDE problems. It essentially boils down to specializing these virtual classes:
  * **Problem**: has knowledge about the discretization used and the system rhs, returns an
  * **Operator**: knows how to relax the solution and can produce a sparse matrix representation of itself
  * **Solver**: given a problem, applies multigrid solver (or an iterative one) to a privately owned solution and is problem & operator agnostic


## Features
  * Extensible design: can solve stencil based discretizations for 1D and 2D problems, given a suitable class implementation
  * Python plotting of experimental results like residual convergence and solve time
  * Convergence tests with problem parameters symbolically generated by sympy
  * Benchmarking of multigrid solver performance with automatic parsing and plotting of results
  * Rich multigrid parameters selection (smoother, cycle, grid transfers...)
  * Jupyter notebook analysis of operator restriction


## Demos
  * `make circuit`: solves a forced RC oscillator circuit and displays the solution
  * `make multilevel`: solves a IsotropicPoisson2D problem with different levels of multigrid and compares their residual converges
  * `make cycle`: compares the convergence of multigrid cycles (full, V, W) on a fixed problem
  * `make convergence_1d`: do a convergence test on the 1D problem defined by the python symbolic generation
  * `make convergence_2d`: do a convergence test on the 2D problem defined by the python symbolic generation
  * `make smoother`: compares the convergence of multigrid solvers on a IsotropicPoisson2D problem varying the smoother and the grid transfer operators

## Benchmarks
  * `make benchmark_omp`: compares multigrid solver on increasingly bigger problems with redblack smoother with different number of openmp threads and plots the speedup
  * `make benchmark_vcycle`: compares non parallelizable gauss-seidel multigrid with parallelizable but weaker red-black and SOR multigrid and plots the speedup
  * `make benchmark_cycle`: compares the performance of same depth V, F and W cycles and plot the speedup
  * `make benchmark_fused`: plots the speedup of kernel fusion optimizations against openmp best multigrid solver
