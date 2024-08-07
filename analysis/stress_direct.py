import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import poisson1D


if __name__ == '__main__':
    x  = sym.symbols('x')
    u  = x**2 * sym.sin(2 * sym.pi * x)

    exact = sym.lambdify(x, u)
    f  = sym.lambdify(x, -sym.diff(u,x,2))

    inf, sup = 0, 1
    boundary = (exact(inf), exact(sup))

    nodes     = 2 ** np.arange(2, 22)
    residuals = np.empty_like(nodes)

    for i, n in enumerate(reversed(nodes)):
        solver = poisson1D.SparseDirectSolver(inf, sup, n, f, boundary)
        solver.solve()

        residuals[i] = np.linalg.norm(solver.residual())

    plt.semilogx(nodes, residuals)
    plt.ylim(0,residuals[-1])
    plt.xlabel('n')
    plt.ylabel('residual norm')
    plt.title("Stress test for direct methods")
    plt.show()
