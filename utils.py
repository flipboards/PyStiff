import numpy as __np
import sympy as __sp
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
