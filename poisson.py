import jax.numpy as jnp
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
        self.xlim = jnp.array(xlim)
        self.size = size

        self.mesh = jnp.linspace(xlim[0], xlim[1], size)

        self.h = jnp.diff(self.xlim) / (size-1)
        self.stencil = Stencil1D(jnp.array([1.0,-2.0,1.0]) / self.h**2)

        self.forcing = forcing
        self.boundary = boundary

        # transfinite interpolation
        x = jnp.linspace(0.0, 1.0, size)
        self.w = boundary(xlim[0]) * (1.0-x) + boundary(xlim[1]) * x + jnp.zeros(size)

        # this rhs solves for the homogeneous dirichlet problem
        self.rhs = self.stencil.apply(self.w, mode='valid') + forcing(self.mesh[1:-1])


    def empty(self):
        return jnp.zeros(self.size-2)


    def assemble(self, homogeneous_solution):
        result = self.w + jnp.pad(homogeneous_solution, pad_width=1)

        return result


class Poisson2D(Problem):
    def __init__(self, xlim: tuple[float, float], ylim: tuple[float, float], size: tuple[int, int], forcing, boundary):
        self.xlim = jnp.array(xlim)
        self.ylim = jnp.array(ylim)
        self.size = jnp.array(size)
        self.inner_size = self.size-2

        self.mesh = jnp.meshgrid(
            jnp.linspace(xlim[0], xlim[1], size[0]),
            jnp.linspace(ylim[0], ylim[1], size[1])
        )

        h = jnp.diff(jnp.stack([self.xlim, self.ylim])).reshape(-1) / (self.size-1)
        self.h = h
        self.stencil = Stencil2D(jnp.array([
            [0.0,                  h[1]**(-2),           0.0       ],
            [h[0]**(-2), -2.0*(h[0]**(-2) + h[1]**(-2)), h[0]**(-2)],
            [0.0,                  h[1]**(-2),           0.0       ]
        ]))

        self.forcing = forcing
        self.boundary = boundary

        # transfinite interpolation
        x = jnp.meshgrid(
            jnp.linspace(0.0, 1.0, size[0]),
            jnp.linspace(0.0, 1.0, size[1])
        )
        g_left   = boundary(xlim[0], self.mesh[1][:,0])
        g_right  = boundary(xlim[1], self.mesh[1][:,0])
        g_top    = boundary(self.mesh[0][0], ylim[0])
        g_bottom = boundary(self.mesh[0][0], ylim[1])

        self.w = (1-x[0])*g_left + x[0]*g_right + (1-x[1])*g_bottom + x[1]*g_top + jnp.zeros(size)

        # this rhs solves for the homogeneous dirichlet problem
        self.rhs = self.stencil.apply(self.w, mode='valid') + forcing(
            self.mesh[0][1:-1, 1:-1],
            self.mesh[1][1:-1, 1:-1]
        )


    def empty(self):
        return jnp.zeros(self.inner_size)


    def assemble(self, homogeneous_solution):
        result = self.w + jnp.pad(homogeneous_solution, pad_width=1)

        return result
