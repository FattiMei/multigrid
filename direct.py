import numpy as np
import scipy


class DirectSolver:
    def __init__(self, problem):
        pass


    def solve(self, rhs):
        pass


# can solve homogeneous Dirichlet problems of every dimensionality, provided a stencil
class SparseSolver(DirectSolver):
    name = 'spsolve'

    def __init__(self, problem):
        self.A = problem.stencil.build_sparse_matrix(problem.size)


    def solve(self, rhs):
        return scipy.sparse.linalg.spsolve(self.A, rhs.flatten()).reshape(rhs.shape)


# the fast poisson solver requires `h`, so I made the constructor take whatever it needs from the problem
class FastSolver(DirectSolver):
    name = 'fast-poisson'

    def __init__(self, problem):
        self.h = problem.h


    # at the moment this works only for 1D Poisson problems, but maybe we could generalize it
    def solve(self, rhs):
        n = len(rhs)

        f_hat = scipy.fft.dst(rhs, type=1)
        eigs = -(2.0 * (1.0 - np.cos(np.pi * np.arange(1,n+1) / (n+1)))) / self.h**2
        u_hat = f_hat / eigs

        return scipy.fft.idst(u_hat, type=1)


class FastSolver2D(FastSolver):
    def solve(self, rhs):
        m,n = rhs.shape
        hx, hy = self.h[0], self.h[1]

        f_hat = scipy.fft.dstn(rhs, type=1)

        eigx = 2.0 * (np.cos(np.pi * np.arange(1,m+1) / (m+1)) - 1.0) / hx**2
        eigy = 2.0 * (np.cos(np.pi * np.arange(1,n+1) / (n+1)) - 1.0) / hy**2

        u_hat = f_hat / (eigx[:,None] + eigy[None,:])

        return scipy.fft.idstn(u_hat, type=1)
