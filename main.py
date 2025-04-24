import numpy as np
import matplotlib.pyplot as plt

from poisson import Poisson1D, Poisson2D
from direct import SparseSolver, FastSolver
from iterative import JacobiSolver, GradientSolver, ConjugateGradient

from time import perf_counter

import sympy as sym
from sympy.abc import x, y


def convergence_test(problem_factory, solvers, nodes, norm=np.linalg.norm):
    err = {solver.name: [] for solver in solvers}

    for n in nodes:
        problem = problem_factory(n)
        exact = problem.boundary(*problem.mesh)

        for solver in solvers:
            homogeneous = solver(problem).solve(problem.rhs)
            sol = problem.assemble(homogeneous)

            err[solver.name].append(
                norm(sol - exact)
            )

    print(err)

    plt.title("$||u - \\tilde{u}||$ - 2nd order expected")
    plt.loglog(nodes, 1.0 / nodes**2, label='$O(N^{-2})$')

    for name, errors in err.items():
        plt.loglog(nodes, errors, label=name)

    plt.legend()
    plt.show()


def iteration_test(problem, solvers, maxit: int, norm=np.linalg.norm):
    assert(maxit > 0)

    plt.title("Convergence behaviour for iterative methods")

    for solver in solvers:
        residuals = np.empty(maxit)
        solver = solver(problem)

        start_time = perf_counter()

        for i in range(maxit):
            residuals[i] = norm(solver.residual())
            solver.step()

        end_time = perf_counter()
        print(solver.name, (end_time - start_time) / maxit)

        plt.semilogy(residuals, label=solver.name)

    plt.legend()
    plt.show()


def discretization_test(u, nodes, norm=np.linalg.norm):
    laplacian = sym.lambdify(
        (x,y),
        sym.diff(u, (x,2)) + sym.diff(u, (y,2))
    )
    u = sym.lambdify((x,y), u)

    err = []

    for n in nodes:
        problem = Poisson2D(
            (0.0, 1.0),
            (0.0, 1.0),
            (n,n),
            lambda x,y: 0.0 * x,
            lambda x,y: 0.0 * x
        )

        stencil = problem.stencil
        exact = laplacian(problem.mesh[0][1:-1,1:-1], problem.mesh[1][1:-1,1:-1])

        err.append(
            norm(exact - stencil.apply(u(*problem.mesh), mode='valid'))
        )

    plt.title("Discretization error")
    plt.loglog(nodes, 1 / nodes**2, label='$O(N^{-2})$')
    plt.loglog(nodes, err, label='err')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    u = sym.sin(x + y) + sym.exp(-x*x + y)
    f = sym.diff(u, (x,2)) + sym.diff(u, (y,2))

    forcing = sym.lambdify((x,y), f)
    boundary = sym.lambdify((x,y), u)

    problem_factory = lambda n: Poisson2D(
        (-1.0, 1.0),
        (-1.0, 1.0),
        (n,n),
        forcing, boundary
    )

    def funky_norm(a):
        b = np.copy(a)
        b[1:-1,1:-1] = 0
        return np.max(np.abs(b))

    convergence_test(
        problem_factory,
        [SparseSolver],
        2 ** np.arange(3,7),
        norm=lambda x: np.max(np.abs(x))
    )

    # iteration_test(
    #     problem_factory(128),
    #     [JacobiSolver, GradientSolver, ConjugateGradient],
    #     maxit=1000
    # )
