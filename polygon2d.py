""" Utils of 2d polygon operations
"""
import numpy as np
import matplotlib.path as mplPath



def _origin_2d(node_coords):
    """ Return the origin point of an 2d polygon. Default is left down.
    """
    return np.argmax(node_coords.dot(np.array([[-1], [-1]])), axis=0)


def _area_2d(node_coords):
    """ Returns the area of an 2d polygon.
    """

    return np.abs(
        np.dot(node_coords[:,0], np.roll(node_coords[:,1], 1)) - 
        np.dot(node_coords[:,1], np.roll(node_coords[:,0], 1))
        ) * 0.5


def _within_2d(point, node_coords):
    """ Returns is a point in a 2d polygon
    """
    boundary = mplPath.Path(node_coords)
    return boundary.contains_point(point)