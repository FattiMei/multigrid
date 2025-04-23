import numpy as np


class IterativeSolver:
    def __init__(self, problem, init=None):
        assert(problem.n > 0)

        if init is None:
            self.x = problem.stencil.apply(np.zeros(problem.n))
        else:
            assert(init.ndim == 1)
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

        center = len(self.stencil.data) // 2
        self.diag = self.stencil.data[center]


    def step(self):
        self.x = self.x + self.residual() / self.diag


class GradientSolver(IterativeSolver):
    name = 'gradient'

    def __init__(self, problem, init=None):
        super().__init__(problem, init)


    def step(self):
        # I'm happy with non optimal implementations, this is just to have another method to test
        residual = self.residual()
        gamma = np.dot(residual, residual) / np.dot(residual, self.stencil.apply(residual, mode='same'))

        self.x = self.x + gamma * residual


class ConjugateGradient(IterativeSolver):
    name = 'cg'

    def __init__(self, problem, init=None):
        super().__init__(problem, init)

        self.p = self.residual()


    def step(self):
        r = self.residual()
        alpha = np.dot(r,r) / np.dot(self.p, self.stencil.apply(self.p, mode='same'))

        self.x = self.x + alpha * self.p
        r_ = r - alpha * self.stencil.apply(self.p, mode='same')

        beta = np.dot(r_,r_) / np.dot(r,r)
        self.p = r_ + beta * self.p
