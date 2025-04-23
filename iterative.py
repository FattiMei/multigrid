import numpy as np
import jax.numpy as jnp


def dot(a, b):
    return jnp.sum(a * b)


class IterativeSolver:
    def __init__(self, problem, init=None):
        if init is None:
            self.x = problem.empty()
        else:
            self.x = init

        self.stencil = problem.stencil
        self.b = problem.rhs


    def step(self):
        pass


    def residual(self):
        return self.b - self.stencil.apply(self.x, mode='same')


class JacobiSolver(IterativeSolver):
    name = 'jacobi'

    def __init__(self, problem, init=None):
        super().__init__(problem, init)
        self.diag = self.stencil.center()


    def step(self):
        self.x = self.x + self.residual() / self.diag


class GradientSolver(IterativeSolver):
    name = 'gradient'

    def __init__(self, problem, init=None):
        super().__init__(problem, init)


    def step(self):
        # I'm happy with non optimal implementations, this is just to have another method to test
        residual = self.residual()
        gamma = dot(residual, residual) / dot(residual, self.stencil.apply(residual, mode='same'))

        self.x = self.x + gamma * residual


class ConjugateGradient(IterativeSolver):
    name = 'cg'

    def __init__(self, problem, init=None):
        super().__init__(problem, init)

        self.p = self.residual()
        self.r = self.p


    def step(self):
        r = self.r
        alpha = dot(r,r) / dot(self.p, self.stencil.apply(self.p, mode='same'))

        self.x = self.x + alpha * self.p

        # to protect from numerical instability we recompute the residual
        r_ = self.residual()

        beta = dot(r_,r_) / dot(r,r)
        self.p = r_ + beta * self.p
        self.r = r_
