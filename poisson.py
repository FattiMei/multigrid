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


def lifting_function(g, X, Y):
    ''' I can't make this work, this lifting makes my solver diverge
    g_x0 = g(0, Y[0, :])      # Left edge
    g_x1 = g(1, Y[-1, :])     # Right edge
    g_y0 = g(X[:, 0], 0)      # Bottom edge
    g_y1 = g(X[:, -1], 1)     # Top edge

    # Corners
    g_00 = g(0, 0)
    g_10 = g(1, 0)
    g_01 = g(0, 1)
    g_11 = g(1, 1)

    # Interpolate across domain
    u_p = (
        (1 - Y) * g_y0[:, np.newaxis] + Y * g_y1[:, np.newaxis] +
        (1 - X) * g_x0[np.newaxis, :] + X * g_x1[np.newaxis, :] -
        ((1 - X) * (1 - Y) * g_00 +
        X * (1 - Y) * g_10 +
        (1 - X) * Y * g_01 +
        X * Y * g_11)
    )
    '''

    return g(X,Y)


class Poisson2D(Problem):
    def __init__(self, xlim: tuple[float, float], ylim: tuple[float, float], size: tuple[int, int], forcing, boundary):
        self.xlim = np.array(xlim)
        self.ylim = np.array(ylim)
        self.size = np.array(size)
        self.inner_size = self.size-2

        x = np.linspace(xlim[0], xlim[1], size[0])
        y = np.linspace(ylim[0], ylim[1], size[1])

        self.mesh = np.meshgrid(x,y, indexing='ij')

        h = np.diff(np.stack([self.xlim, self.ylim])).reshape(-1) / (self.size-1)
        self.h = h
        self.stencil = Stencil2D(np.array([
            [0.0,                  h[1]**(-2),           0.0       ],
            [h[0]**(-2), -2.0*(h[0]**(-2) + h[1]**(-2)), h[0]**(-2)],
            [0.0,                  h[1]**(-2),           0.0       ]
        ]))

        self.forcing = forcing
        self.boundary = boundary

        # transfinite interpolation, Coons patch doesn't work
        self.w = lifting_function(boundary, self.mesh[0], self.mesh[1])

        # this rhs solves for the homogeneous dirichlet problem
        self.rhs = forcing(
            self.mesh[0][1:-1, 1:-1],
            self.mesh[1][1:-1, 1:-1]
        ) - self.stencil.apply(self.w, mode='valid')


    def empty(self):
        return np.zeros(self.inner_size)


    def assemble(self, homogeneous_solution):
        return self.w + np.pad(homogeneous_solution, pad_width=1)
