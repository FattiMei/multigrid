import numpy as np
import sympy as sym
from sympy.abc import x
import matplotlib.pyplot as plt

from poisson import Poisson1D
from direct import SparseSolver, FastSolver


def convergence_test(method, nodes, norm=np.linalg.norm):
    u = sym.sin(2*x) * 10.0 * sym.exp(-x*x) + x**2
    f = sym.diff(u, (x,2))

    inf, sup = 0.0, 10.0
    forcing = sym.lambdify(x, f)
    boundary = sym.lambdify(x, u)

    err = []

    for n in nodes:
        problem = Poisson1D(n, inf, sup, forcing, boundary)
        stencil = problem.stencil
        rhs = problem.rhs

        homogeneous = method(problem).solve(rhs)
        sol = problem.assemble(homogeneous)

        err.append(
            norm(sol - boundary(problem.mesh))
        )

    return err


if __name__ == '__main__':
    nodes = 10 ** np.arange(1, 6)

    plt.loglog(nodes, 1.0 / nodes**2, label='$O(N^{-2})$')

    for solver in [SparseSolver, FastSolver]:
        err = convergence_test(solver, nodes, norm=lambda x: np.max(np.abs(x)))
        plt.loglog(nodes, err, label=solver.name)

    plt.legend()
    plt.show()
