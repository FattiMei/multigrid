import numpy as np
import scipy


class Stencil1D:
    ndim = 1

    def __init__(self, stencil: np.ndarray):
        assert(stencil.ndim == 1)
        assert(stencil.size > 0)

        self.data = stencil


    # variation of the stencil class will modify the convolution implementation
    def apply(self, x, mode='valid'):
        return np.convolve(x, self.data, mode)


    def build_sparse_matrix(self, n: int):
        assert(self.data.size == 3)
        assert(n > 2)

        return scipy.sparse.diags(
            self.data,
            [-1,0,1],
            shape=(n-2,n-2)
        ).tocsr()


class Stencil2D:
    ndim = 2

    def __init__(self, stencil: np.ndarray):
        assert(stencil.ndim == 2)
        assert(np.all(stencil.shape > 0))

        self.data = stencil


    def apply(self, x, mode='valid'):
        pass


    def build_sparse_matrix(self, rows: int, cols: int):
        assert(self.data.shape == (3,3))
        assert(rows > 0)
        assert(cols > 0)

        pass
