import numpy as __np

def xyz2array(xyz_data):
    """ Convert x,y,z,... -> np.array(x,y,z)
    """
    return __np.array(xyz_data.split(',')).astype(__np.float)


