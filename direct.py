import numpy as np
import scipy


class DirectSolver:
    def __init__(self, n: int, stencil):
        pass


    def solve(self, rhs):
        pass


class SparseSolver(DirectSolver):
    def __init__(self, n: int, stencil):
        self.name = 'sparse'

        # Assumes a 3 point stencil and stores only internal points
        self.A = scipy.sparse.diags(stencil, [-1,0,1], shape=(n-2,n-2)).tocsr()


    def solve(self, rhs):
        return scipy.sparse.linalg.spsolve(self.A, rhs)


class FastSolver(DirectSolver):
    # REMARK: works only for the homogeneous regular 1D Poisson problem
    def __init__(self, n: int, stencil):
        self.name = 'fast-poisson'


    def solve(self, rhs):
        pass
