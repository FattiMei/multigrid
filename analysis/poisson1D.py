import numpy as np
import sympy as sym
import scipy
import unittest


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


    def solve(self):
        pass


class MatrixFreeSolver(Poisson1D):
    def __init__(self, inf, sup, n, f, boundary):
        super().__init__(inf, sup, n, f, boundary)
        self.it = 0


    def action(self, x):
        y = np.empty_like(x)
        y[0] = x[0]
        y[-1] = x[-1]

        for i in range(1,self.n-1):
            y[i] = (2*x[i] - x[i-1] - x[i+1]) / self.h**2

        return y


    def residual(self):
        return self.b - self.action(self.u)


    def step(self):
        self.it = self.it + 1


    def solve(self, tol, maxiter):
        self.status = 'maxiter'
        self.r = self.residual()

        for _ in range(maxiter):
            self.step()
            self.r = self.residual()

            if np.linalg.norm(self.r) < tol:
                self.status = 'ok'
                break


class JacobiSolver(MatrixFreeSolver):
    def __init__(self, inf, sup, n, f, boundary):
        super().__init__(inf, sup, n, f, boundary)
        self.label = 'jacobi'


    def step(self):
        super().step()

        y = np.zeros_like(self.u)
        y[0] = self.u[0]
        y[-1] = self.u[-1]

        for i in range(1,self.n-1):
            y[i] = (self.h**2 * self.b[i] + self.u[i-1] + self.u[i+1]) / 2

        self.u = y


class GaussSeidelSolver(MatrixFreeSolver):
    def __init__(self, inf, sup, n, f, boundary):
        super().__init__(inf, sup, n, f, boundary)
        self.label = 'gseidel'


    def step(self):
        super().step()

        for i in range(1,self.n-1):
            self.u[i] = (self.h**2 * self.b[i] + self.u[i-1] + self.u[i+1]) / 2


class CgSolver(MatrixFreeSolver):
    def __init__(self, inf, sup, n, f, boundary):
        super().__init__(inf, sup, n, f, boundary)
        self.label = 'cg'
        self.r = self.residual()
        self.d = self.r


    def step(self):
        super().step()

        y      = self.action(self.d)

        alpha  = np.dot(self.r, self.r) / np.dot(self.d, y)
        old    = self.r

        self.u = self.u + alpha * self.d
        self.r = self.residual()

        beta   = np.dot(self.r, self.r) / np.dot(old,old)
        self.d = self.r + beta * self.d


    # to remove the unnecessary recomputing of the residual, there is this little duplication in the solving logic
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


    def solve(self, *args):
        self.u      = np.linalg.solve(self.A, self.b)
        self.status = 'ok'


class SparseDirectSolver(Poisson1D):
    def __init__(self, inf, sup, n, f, boundary):
        super().__init__(inf, sup, n, f, boundary)
        self.label = 'sparse'

        # I need to build the matrix A as sparse because memory is a problem
        h        = self.h
        diagonal = np.concatenate([np.ones(1), 2*np.ones(n-2)/h**2, np.ones(1)])
        upper    = np.concatenate([np.zeros(1), -np.ones(n-2)/h**2])
        lower    = np.concatenate([-np.ones(n-2)/h**2, np.zeros(1)])

        self.A   = scipy.sparse.diags([diagonal, upper, lower], [0, 1, -1], format='csr')


    def dense_repr(self):
        return self.A.toarray()


    def residual(self):
        return self.b - self.A.dot(self.u)


    def solve(self):
        self.u = scipy.sparse.linalg.spsolve(self.A, self.b)
        self.status = 'ok'


class TestPoissonSolvers(unittest.TestCase):
    def setUp(self):
        x  = sym.symbols('x')
        u  = x**2 * sym.sin(2 * sym.pi * x)

        exact = sym.lambdify(x, u)
        f  = sym.lambdify(x, -sym.diff(u,x,2))

        n        = 100
        inf, sup = 0, 1
        boundary = (exact(inf), exact(sup))

        self.dense  = DirectSolver      (inf, sup, n, f, boundary)
        self.sparse = SparseDirectSolver(inf, sup, n, f, boundary)
        self.mfree  = MatrixFreeSolver  (inf, sup, n, f, boundary)


    def test_sparse_repr(self):
        self.assertTrue(np.allclose(self.dense.A, self.sparse.dense_repr()))


    def test_direct_solvers(self):
        self.dense.solve()
        self.sparse.solve()

        self.assertTrue(np.allclose(self.dense.u, self.sparse.u))


    def test_matrix_free_action(self):
        x = np.random.rand(self.dense.u.shape[0])

        self.assertTrue(np.allclose(self.dense.A @ x, self.sparse.A.dot(x)))
        self.assertTrue(np.allclose(self.dense.A @ x, self.mfree.action(x)))


    def test_residual_calculation(self):
        x = np.random.rand(self.dense.u.shape[0])

        self.dense.u = x
        self.mfree.u = x

        self.assertTrue(np.allclose(self.dense.residual(), self.mfree.residual()))


if __name__ == '__main__':
    unittest.main()
