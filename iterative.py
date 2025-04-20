import numpy as np


class IterativeSolver:
    def __init__(self, n: int, stencil, init=None):
        assert(n > 0)

        self.n = n
        self.stencil = stencil

        if init is None:
            self.x = np.zeros(n)
        else:
            assert(init.ndim == 1)
            assert(init.shape[0] == n)

            self.x = init


    def step():
        pass


# TODO: implement all the easy black box solvers
class JacobiSolver(IterativeSolver):
    def __init__(self, n: int, stencil, init=None):
        super().__init__()
        self.name = 'jacobi'


    def step():
        pass
