import numpy as np
from datetime import datetime
from dateutil import parser

from specklepy.exceptions import SpecklepyTypeError


def weighted_mean(data, axis=None, weights=None, vars=None):
    """Combine data frames considering weights and variances.

    Args:
        data (np.ndarray):
        weights (np.ndarray, optional):
        vars (np.ndarray, optional):

    Returns:
        weighted_mean (np.ndarray):
        err (np.ndarray, optional):
            Returned only if vars is provided.
    """

    if vars is not None:
        vars = np.ma.masked_values(vars, 0.0)
        if weights is None:
            weights = np.divide(1, vars)
        else:
            weights = np.divide(weights, vars)

    mean, weights = np.average(data, axis=axis, weights=weights, returned=True)

    weights = np.ma.masked_values(weights, 0.0)  # Avoid division by zero

    return mean, np.divide(1, weights)


def time_difference(t0, times):
    """Compute time difference(s) from one time to others.

    Args:
        t0 (str, datetime):
            Time stamp whose difference to times is computed.
        times (str, datetime, list):
            Times to which the time difference is computed.

    Returns:
        timedeltas (float, np.array):
            Time difference(s) between from t0 to times (in seconds).
    """

    # Type check
    if isinstance(t0, str):
        t0 = parser.parse(t0)
    elif isinstance(t0, datetime):
        pass
    else:
        raise SpecklepyTypeError('time_difference', 't0', type(t0), 'str or datetime')

    if isinstance(times, str):
        times = parser.parse(times)
    elif isinstance(times, (datetime, list, np.ndarray)):
        pass
    else:
        raise SpecklepyTypeError('time_difference', 'times', type(times), 'str, datetime, or list')

    # Compute time deltas in units of seconds
    if isinstance(times, (list, np.ndarray)):
        timedeltas = np.empty((len(times)))
        for i, o in enumerate(times):
            if isinstance(o, str):
                o = parser.parse(o)
            timedeltas[i] = (o - t0).total_seconds()
        return timedeltas
    else:
        return (times - t0).total_seconds()


def get_distance_weights(distances, norm=True):
    """

    Args:
        distances (np.ndarray):
            Spatial or time distances.
        norm (bool):
            Normalize to unity.

    Returns:
        weights (np.ndarray):
            Weights, normed if requested.
    """

    weights = np.square(np.divide(1, distances))

    if norm:
        weights = np.divide(weights, np.sum(weights))

    return weights


def combine_masks(*arrays):

    def get_mask(arr):
        if arr is None:
            return None
        try:
            return arr.mask
        except AttributeError:
            return np.zeros(arr.shape, dtype=bool)

    mask = get_mask(arrays[0])
    try:
        for array in arrays[1:]:
            mask = np.logical_or(mask, get_mask(array))
        return mask
    except IndexError:
        return mask
