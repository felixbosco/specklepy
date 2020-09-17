import numpy as np


def rad2deg(val):
    """Transform an angle from units of radians into degrees."""
    return np.rad2deg(val)


def rad2mas(val):
    """Transform an angle from units of radians into milli-arcseconds."""
    return np.rad2deg(val) * 60**2 * 1000
