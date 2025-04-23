import os
import scipy

if os.environ.get('USE_JAX') is None:
    import numpy as np
    from scipy.signal import convolve2d
else:
    import jax.numpy as np
    from jax.scipy.signal import convolve2d


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
        assert(self.data.size == 3)
        assert(n > 2)

        # idx  = np.arange(n)
        # rows = []
        # cols = []
        # data = []
        #
        # for i in range(self.data.shape[0]):
        #     offset = self.data.shape[0]//2 - i
        #     neighbours = idx - offset
        #
        #     valid = np.argwhere((neighbours >= 0) & (neighbours < n)).flatten()
        #
        #     rows.append(idx[valid])
        #     cols.append(neighbours[valid])
        #     data.append(self.data[i] * np.ones(len(valid)))
        #
        # return scipy.sparse.coo_array(
        #     (np.concat(data), (np.concat(rows), np.concat(cols)))
        # ).tocsr()

        return scipy.sparse.diags(
            self.data,
            [-1,0,1],
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


    def center(self):
        return self.data[self.data.shape[0] // 2, self.data.shape[1] // 2]
