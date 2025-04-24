import os
import scipy
import numpy as np

from scipy.signal import convolve2d


class Stencil:
    def __init__(self):
        pass


    def apply(self, x, mode):
        pass


    def build_sparse_matrix(self, size):
        pass


    def center(self):
        pass


class Stencil1D(Stencil):
    ndim = 1

    def __init__(self, stencil: np.ndarray):
        assert(stencil.ndim == 1)
        assert(stencil.size > 0)

        self.data = stencil


    # variation of the stencil class will modify the convolution implementation
    def apply(self, x, mode='valid'):
        return np.convolve(x, self.data, mode)


    def build_sparse_matrix(self, n: int):
        size = len(self.data)
        offsets = np.arange(size) - size//2

        return scipy.sparse.diags(
            self.data,
            offsets,
            shape=(n-2,n-2)
        ).tocsr()


    def center(self):
        return self.data[self.data.shape[0] // 2]


class Stencil2D(Stencil):
    ndim = 2

    def __init__(self, stencil: np.ndarray):
        assert(stencil.ndim == 2)
        assert(np.all(np.array(stencil.shape) > 0))

        self.data = stencil


    def apply(self, x, mode='valid'):
        return convolve2d(x, self.data, mode)


    def build_sparse_matrix(self, size: tuple[int, int]):
        assert(self.data.shape == (3,3))

        n,m = size
        N = (n-2) * (m-2)

        offsets = np.array([-1,0,1])

        return scipy.sparse.diags(
            self.data.flatten(),
            (offsets + n*offsets[:, None]).flatten(),
            shape=(N,N)
        ).tocsr()


    def center(self):
        return self.data[self.data.shape[0] // 2, self.data.shape[1] // 2]
