import numpy as np
import matplotlib.pyplot as plt

from poisson import Poisson1D
from direct import SparseSolver, FastSolver


def convergence_test(problem_factory, solvers, nodes, norm=np.linalg.norm):
    err = {solver.name: [] for solver in solvers}

    for n in nodes:
        problem = problem_factory(n)
        exact = problem.boundary(problem.mesh)

        for solver in solvers:
            homogeneous = solver(problem).solve(problem.rhs)
            sol = problem.assemble(homogeneous)

            err[solver.name].append(
                norm(sol - exact)
            )

    plt.title("Convergence test - 2nd order expected")
    plt.loglog(nodes, 1.0 / nodes**2, label='$O(N^{-2})$')

    for name, errors in err.items():
        plt.loglog(nodes, errors, label=name)

    plt.legend()
    plt.show()


if __name__ == '__main__':
    import sympy as sym
    from sympy.abc import x

    u = sym.sin(2*x) * 10.0 * sym.exp(-x*x) + x**2
    f = sym.diff(u, (x,2))

    inf, sup = 0.0, 10.0
    forcing = sym.lambdify(x, f)
    boundary = sym.lambdify(x, u)

    problem_factory = lambda n: Poisson1D(n, inf, sup, forcing, boundary)
    nodes = 10 ** np.arange(1,6)

    convergence_test(
        lambda n: Poisson1D(n, inf, sup, forcing, boundary),
        [SparseSolver, FastSolver],
        nodes
    )
