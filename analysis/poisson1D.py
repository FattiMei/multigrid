import numpy as np


class Poisson1D:
    def __init__(self, inf, sup, n, f, boundary):
        self.mesh  = np.linspace(inf, sup, n)
        self.n     = n
        self.h     = (sup - inf) / (n-1)

        self.u     = np.zeros(n)
        self.u[0]  = boundary[0]
        self.u[-1] = boundary[1]

        self.b     = f(self.mesh)
        self.b[0]  = self.u[0]
        self.b[-1] = self.u[-1]


    def residual(self):
        pass


    def step(self):
        pass


    def solve(self, tol, maxiter):
        self.status = 'maxiter'

        for _ in range(maxiter):
            self.step()

            if np.linalg.norm(self.r) < tol:
                self.status = 'ok'
                break


class DirectSolver(Poisson1D):
    def __init__(self, inf, sup, n, f, boundary):
        super().__init__(inf, sup, n, f, boundary)
        self.label = 'direct'

        A = np.zeros((n,n))
        for i in range(1,n-1):
            A[i,i-1] = -1
            A[i,i]   =  2
            A[i,i+1] = -1

        A = A / self.h ** 2
        A[0,0] = 1
        A[n-1,n-1] = 1

        self.A = A


    def residual(self):
        return self.b - self.A @ self.u


    def solve(self, tol, maxiter):
        self.u      = np.linalg.solve(self.A, self.b)
        self.status = 'ok'
