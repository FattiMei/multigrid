import numpy as np
import scipy


class DirectSolver:
    def __init__(self, problem):
        pass


    def solve(self, rhs):
        pass


class SparseSolver(DirectSolver):
    name = 'spsolve'

    def __init__(self, problem):

        # Assumes a 3 point stencil and stores only internal points
        self.A = scipy.sparse.diags(
            problem.stencil,
            [-1,0,1],
            shape=(problem.n-2,problem.n-2)
        ).tocsr()


    def solve(self, rhs):
        return scipy.sparse.linalg.spsolve(self.A, rhs)


# the fast poisson solver requires `h`, so I made the constructor take whatever it needs from the problem
class FastSolver(DirectSolver):
    name = 'fast-poisson'

    def __init__(self, problem):
        self.h = problem.h


    def solve(self, rhs):
        n = len(rhs)

        f_hat = scipy.fft.dst(rhs, type=1)
        eigs = -(2.0 * (1.0 - np.cos(np.pi * np.arange(1,n+1) / (n+1)))) / self.h**2
        u_hat = f_hat / eigs

        return scipy.fft.idst(u_hat, type=1)
