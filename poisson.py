import numpy as np
from stencil import Stencil1D


class Poisson1D:
    def __init__(self, n: int, inf: float, sup: float, forcing, boundary):
        assert(n > 3)
        assert(inf < sup)

        self.n = n
        self.h = (sup-inf) / (n-1.0)
        self.mesh = np.linspace(inf, sup, n)
        self.stencil = Stencil1D(np.array([1.0,-2.0,1.0]) / self.h**2)

        self.forcing = forcing
        self.boundary = boundary

        # transfinite interpolation
        x = np.linspace(0.0, 1.0, n)
        self.w = boundary(inf) * (1.0-x) + boundary(sup) * x

        # this rhs solves for the homogeneous dirichlet problem
        self.rhs = forcing(self.mesh[1:-1]) + self.stencil.apply(self.w, mode='valid')


    def assemble(self, homogeneous_solution):
        # we need a copy, otherwise multiple solvers will corrupt this piece of data
        result = np.copy(self.w)
        result[1:-1] += homogeneous_solution

        return result
