import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import poisson1D


MAXITER = 1000
TOL     = 1e-6


if __name__ == '__main__':
    x  = sym.symbols('x')
    u  = x**2 * sym.sin(2 * sym.pi * x)

    exact = sym.lambdify(x, u)
    f  = sym.lambdify(x, -sym.diff(u,x,2))

    inf, sup = 0, 1
    boundary = (exact(inf), exact(sup))

    nodes = [2**i for i in range(2,8)]
    solver_classes = [
        poisson1D.JacobiSolver,
        poisson1D.GaussSeidelSolver,
        poisson1D.CgSolver,
    ]

    iterations = {
        'jacobi' : [],
        'gseidel': [],
        'cg'     : [],
    }

    for n in nodes:
        for solver in solver_classes:
            s = solver(inf, sup, n, f, boundary)
            s.solve(TOL, MAXITER)
            iterations[s.label].append(s.it)

    for label, it in iterations.items():
        plt.plot(nodes, it, label=label)

    plt.title(f"Iterations until norm(r) < {TOL}")
    plt.ylim(0,MAXITER-1)
    plt.xlabel('problem size')
    plt.ylabel('it')
    plt.legend()
    plt.show()
