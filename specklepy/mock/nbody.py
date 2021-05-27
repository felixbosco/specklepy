import numpy as np
import warnings

from astropy.constants import G
from astropy.units import Unit, Quantity

from specklepy.logging import logger


def vector_lengths(vectors, axis=-1):
    return np.sqrt(np.sum(np.square(vectors), axis=axis))


def acceleration(positions, masses, position_unit='pc', mass_unit='M_sun'):
    """Compute vectors of the gravitational acceleration by the other bodies.

    Args:
        positions (np.ndarray):
            Array of position vectors. Must have shape `(n, 3)`.
        masses (np.ndarray):
            Array of masses of the bodies. Must have shape `(n,)`.
        position_unit (str, optional):
            Unit of the `positions` array. Default is `'pc'`.
        mass_unit (str, optional):
            Unit of the `masses` array. Default is `'M_sun'`, i.e. solar masses.

    Returns:
        accelerations (np.ndarray):
            Array of acceleration vectors in units of m s-2. Has the same shape as `positions`.
    """

    # Derive gravitational constant, such that accelerations is in m s-2
    gravitational_constant = (G * Unit(f"{position_unit}-2 {mass_unit}")).decompose().value

    # Derive distance vectors and normalize
    delta_positions = positions[np.newaxis, :, :] - positions[:, np.newaxis, :]
    absolute_radius = vector_lengths(delta_positions)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        normed = np.divide(delta_positions, np.power(absolute_radius, 3)[:, :, np.newaxis])
    del delta_positions
    del absolute_radius

    # Compute acceleration vectors
    accelerations = gravitational_constant * np.nansum(np.multiply(masses[np.newaxis, :, np.newaxis], normed), axis=1)

    return accelerations


def orbital_velocity(positions, masses, position_unit='pc', mass_unit='M_sun', rotation_axis=None):
    """Compute orbital velocity vectors that match roughly the gravitational acceleration by the other bodies.

    Args:
        positions (np.ndarray):
            Array of position vectors. Must have shape `(n, 3)`.
        masses (np.ndarray):
            Array of masses of the bodies. Must have shape `(n,)`.
        position_unit (str, optional):
            Unit of the `positions` array. Default is `'pc'`.
        mass_unit (str, optional):
            Unit of the `masses` array. Default is `'M_sun'`, i.e. solar masses.
        rotation_axis (np.ndarray):
            Unit vector of the rotation axis. Falls back to to z-axis, i.e. `np.array([0, 0, 1])`, if not provided.

    Returns:
        velocity_vectors (np.ndarray):
            Array of velocity vectors in units of m s-1. Has the same shape as `positions`.
    """

    if rotation_axis is None:
        logger.info(f"Setting rotation axis to default: z-axis")
        rotation_axis = np.array([0, 0, 1])

    accelerations = acceleration(positions, masses, position_unit=position_unit, mass_unit=mass_unit)
    abs_accelerations = vector_lengths(accelerations)
    abs_radii = vector_lengths(positions) * Quantity(f"1 {position_unit}").to('m').value
    scales = np.sqrt(np.divide(abs_radii, abs_accelerations))

    # Compute the velocity vectors
    velocity_vectors = np.cross(rotation_axis, scales[:, np.newaxis] * accelerations)

    return velocity_vectors
