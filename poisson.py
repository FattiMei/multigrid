import numpy as np
from stencil import Stencil1D, Stencil2D


class Problem:
    def __init__(self):
        pass


    def empty(self):
        pass


    def assemble(self, homogeneous_solution):
        pass


class Poisson1D(Problem):
    def __init__(self, xlim: tuple[float, float], size: tuple[int], forcing, boundary):
        self.xlim = np.array(xlim)
        self.size = np.array(size)

        self.mesh = np.linspace(xlim[0], xlim[1], size)

        self.h = np.diff(self.xlim) / (size-1)
        self.stencil = Stencil1D(np.array([1.0,-2.0,1.0]) / self.h**2)

        self.forcing = forcing
        self.boundary = boundary

        # transfinite interpolation
        x = np.linspace(0.0, 1.0, size)
        self.w = boundary(xlim[0]) * (1.0-x) + boundary(xlim[1]) * x + np.zeros(size)

        # this rhs solves for the homogeneous dirichlet problem
        self.rhs = forcing(self.mesh[1:-1]) - self.stencil.apply(self.w, mode='valid')


    def empty(self):
        return np.zeros(self.size-2)


    def assemble(self, homogeneous_solution):
        result = self.w + np.pad(homogeneous_solution, pad_width=1)

        assert(result[0] == self.w[0])
        assert(result[-1] == self.w[-1])

        return result


class Poisson2D(Problem):
    def __init__(self, xlim: tuple[float, float], ylim: tuple[float, float], size: tuple[int, int], forcing, boundary):
        self.xlim = np.array(xlim)
        self.ylim = np.array(ylim)
        self.size = np.array(size)
        self.inner_size = self.size-2

        self.mesh = np.meshgrid(
            np.linspace(xlim[0], xlim[1], size[0]),
            np.linspace(ylim[0], ylim[1], size[1])
        )

        h = np.diff(np.stack([self.xlim, self.ylim])).reshape(-1) / (self.size-1)
        self.h = h
        self.stencil = Stencil2D(np.array([
            [0.0,                  h[1]**(-2),           0.0       ],
            [h[0]**(-2), -2.0*(h[0]**(-2) + h[1]**(-2)), h[0]**(-2)],
            [0.0,                  h[1]**(-2),           0.0       ]
        ]))

        self.forcing = forcing
        self.boundary = boundary

        # transfinite interpolation (Coons patch)
        # u = np.linspace(0,1,size[0])
        # v = np.linspace(0,1,size[1])
        # uu, vv = np.meshgrid(u,v)
        #
        # boundary_left = boundary(xlim[0], np.linspace(ylim[0], ylim[1], size[1]))
        # boundary_right = boundary(xlim[1], np.linspace(ylim[0], ylim[1], size[1]))
        # boundary_top = boundary(np.linspace(xlim[0], xlim[1], size[0]), ylim[1])
        # boundary_bottom = boundary(np.linspace(xlim[0], xlim[1], size[0]), ylim[0])
        #
        # L = (1-uu) * boundary_left + uu * boundary_right
        # B = (1-vv) * boundary_bottom + vv * boundary_top
        #
        # corner_blend = (
        #     (1-uu) * (1-vv) * boundary_bottom[0] +
        #     uu * (1-vv) * boundary_bottom[-1] +
        #     (1-uu) * vv * boundary_top[0] +
        #     uu * vv * boundary_top[-1]
        # )
        # self.w = L + B - corner_blend

        self.w = boundary(*self.mesh)

        # this rhs solves for the homogeneous dirichlet problem
        self.rhs = forcing(
            self.mesh[0][1:-1, 1:-1],
            self.mesh[1][1:-1, 1:-1]
        ) - self.stencil.apply(self.w, mode='valid')


    def empty(self):
        return np.zeros(self.inner_size)


    def assemble(self, homogeneous_solution):
        return self.w + np.pad(homogeneous_solution, pad_width=1)
