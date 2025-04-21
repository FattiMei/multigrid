import numpy as np
import sympy as sym
from sympy.abc import x

from poisson import Poisson1D
from direct import SparseSolver


def convergence_test(method, nodes, norm=np.linalg.norm):
    u = sym.sin(2*x) + sym.exp(-x*x)
    f = sym.diff(u, (x,2))

    inf, sup = 0.0, 1.0
    forcing = sym.lambdify(x, f)
    boundary = sym.lambdify(x, u)

    err = []

    for n in nodes:
        problem = Poisson1D(n, inf, sup, forcing, boundary)
        stencil = problem.stencil
        rhs = problem.rhs

        homogeneous = method(n, stencil).solve(rhs)
        sol = problem.assemble(homogeneous)

        err.append(
            norm(sol - boundary(problem.mesh))
        )

    return err


if __name__ == '__main__':
    nodes = 10 ** np.arange(1, 6)
    err = convergence_test(SparseSolver, nodes)

    plt.loglog(nodes, err, label='spsolve')
    plt.loglog(nodes, 1.0 / nodes**2, label='$O(N^{-2})$')

    plt.legend()
    plt.show()
