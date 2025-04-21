import numpy as np
import scipy
import direct
import iterative


class Poisson1D:
    def __init__(self, n: int, inf: float, sup: float, forcing, boundary):
        assert(n > 3)
        assert(inf < sup)

        h = (sup-inf) / (n-1.0)
        self.mesh = np.linspace(inf, sup, n)
        self.stencil = np.array([1.0,-2.0,1.0]) / (h*h)

        # transfinite interpolation
        x = np.linspace(inf, sup, n)
        self.w = boundary(inf) * (1.0-x) + boundary(sup) * x

        # this rhs solves for the homogeneous dirichlet problem
        self.rhs = forcing(self.mesh[1:-1]) + np.convolve(
            self.w,
            self.stencil,
            mode='valid'
        )


    def assemble(self, homogeneous_solution):
        result = self.w
        result[1:-1] += homogeneous_solution

        return result
