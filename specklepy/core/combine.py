import numpy as np

def weighted_combine(data, axis=None, weights=None, vars=None):
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
        if weights is None:
            weights = np.divide(1, vars)
        else:
            weights = np.divide(weights, vars)

    mean, weights = np.average(data, axis=axis, weights=weights, returned=True)

    weights = np.ma.masked_values(weights, 0.0) # Avoid division by zero

    return mean, np.divide(1, weights)
