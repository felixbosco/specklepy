from IPython import embed
import numpy as np

from astropy import stats


def sigma_clip(cube, **kwargs):
    """

    Args:
        cube:
        **kwargs:

    Returns:

    """
    out = np.ma.empty_like(cube)

    for row_index in range(cube.shape[1]):
        out[:, row_index, :] = stats.sigma_clip(cube[:, row_index, :], axis=0, masked=True, **kwargs)

    return out


def sigma_clipped_stats(cube, **kwargs):
    """

    Args:
        cube:
        **kwargs:

    Returns:

    """
    number_rows = cube.shape[1]
    frame_shape = cube.shape[1:]
    mean = np.empty(frame_shape)
    median = np.empty(frame_shape)
    std = np.empty(frame_shape)

    for row_index in range(number_rows[:2]):
        row_stats = stats.sigma_clipped_stats(cube[:, row_index, :], axis=0, **kwargs)
        mean[:, row_index, :], median[:, row_index, :], std[:, row_index, :] = row_stats

    return mean, median, std
