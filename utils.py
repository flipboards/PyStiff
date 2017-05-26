import numpy as __np
import sympy as __sp
import multiprocessing as __mp
from sympy.abc import x

def xyz2array(xyz_data):
    """ Convert x,y,z,... -> np.array(x,y,z)
    """
    return __np.array(xyz_data.split(',')).astype(__np.float)


def str_mat(matrix, format=False, symbol=False):
    """ format matrix
    """

    if isinstance(matrix, list) or len(matrix.shape) == 1:

        return '\n'.join(['\t% .4e' % a for a in matrix])    
    
    elif len(matrix.shape) == 2:

        s = ''

        if symbol:
            for i in range(matrix.shape[0]):
                s += '\t' + '\t'.join([str(a.evalf(4)) for a in matrix[i,:]]) + '\n'

        elif not format:

            for i in range(matrix.shape[0]):
                s += '\t'.join([str(a) for a in matrix[i,:]]) + '\n'
        else:
            for i in range(matrix.shape[0]):
                s += '\t' + '\t'.join(['% .4e' % a for a in matrix[i,:]]) + '\n'
        return s


def sym2float(sym):
    return __sp.lambdify((x), sym)(0)

def __dot(x, y):
    return __np.dot(x, y)


def mpdot(a, b):

    pool = __mp.Pool(4)

    tmp = []

    for i in range(a.shape[0]):
        tmp.append([])
        for j in range(b.shape[1]):
            tmp[-1].append(pool.apply_async(__dot, args=(a[i, :], b[:, j])))

    pool.close()
    pool.join()

    result = __np.zeros((a.shape[0], b.shape[1]), dtype=object)
    for i in range(a.shape[0]):
        for j in range(b.shape[1]):
            result[i, j] = tmp[i][j].get()

    return result