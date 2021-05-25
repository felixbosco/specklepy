import numpy as np
import warnings

from astropy.constants import G
from astropy.units import Unit


def acceleration(positions, masses, position_unit='pc', mass_unit='M_sun'):
    """

    Args:
        positions (np.ndarray):
            Array of position vectors. Must have shape `(n, 3)`.
        masses (np.ndarray):
            Array of masses of the bodies. Must have shape `(n, )`.
        position_unit (str, optional):
            Unit of the `positions` array. Default is `'pc'`.
        mass_unit (str, optional):
            Unit of the `masses` array. Default is `'M_sun'`, i.e. solar masses.

    Returns:
        accelerations (np.ndarray):
            Array of accelerations in units of m s-2. Has the same shape as `positions`.
    """

    # Derive gravitational constant, such that accelerations is in m s-2
    gravitational_constant = (G * Unit(f"{position_unit}-2 {mass_unit}")).decompose().value

    # Derive distance vectors and normalize
    delta_positions = positions[np.newaxis, :, :] - positions[:, np.newaxis, :]
    norm_squared = np.sum(np.square(delta_positions), axis=-1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        normed = np.divide(delta_positions, np.power(norm_squared, 3/2)[:, :, np.newaxis])
    del delta_positions
    del norm_squared

    # Compute acceleration vectors
    accelerations = gravitational_constant * np.nansum(np.multiply(masses[np.newaxis, :, np.newaxis], normed), axis=1)

    return accelerations
