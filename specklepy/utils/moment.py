import numpy as np


def moment_2d(array, var=None, weights=None):
    """Compute zeroth and first order moments of an array.

    Uncertainties are `None` if no variance is provided.

    Arguments:
        array (array-like):
            .
        var (array-like):
            Array of variances with the same shape as `array` argument. If not provided, all uncertainties will be
            `None`-type.
        weights (array-like):
            Weights to apply during the computation. If not provided, all entries of the `array` and `var` input are
            weighted equally.

    Returns:
        moment0 (float):
            Zeroth order Moment.
        std_moment0 (float):
            Uncertainty of `moment0`.
        moment1_x (float):
            First order moment in x-direction.
        std_moment1_x (float):
            Uncertainty of `moment1_x`. Will be None if `var` is not provided.
        moment1_y (float):
            First order moment in y-direction.
        std_moment1_y (float):
            Uncertainty of `moment1_y`. Will be None if `var` is not provided.
    """

    # Apply weights to array and variance map
    try:
        weighted_array = np.multiply(weights, array)
    except TypeError:
        weighted_array = array
    try:
        weighted_var = np.multiply(weights, var)
    except TypeError:
        weighted_var = var

    # Initialize arrays
    y, x = np.mgrid[:array.shape[0], :array.shape[1]]

    # Compute zeroth moment
    moment0 = np.sum(weighted_array)
    try:
        var_moment0 = np.sum(np.square(weighted_var))
        std_moment0 = np.sqrt(var_moment0)
    except TypeError:
        std_moment0 = None

    # Compute first moments
    try:
        moment1_x = np.average(x, weights=weighted_array, axis=(0, 1))
        moment1_y = np.average(y, weights=weighted_array, axis=(0, 1))
    except ZeroDivisionError as e:
        if (np.array(array.shape) == 0).any():
            raise ValueError(f"Input array with shape {array.shape} has no extent and cannot be averaged!")
        else:
            raise e

    # Compute uncertainties
    try:
        var_moment1_x = np.sum(np.square(np.multiply(weighted_var, np.subtract(x, moment1_x)))) / np.square(moment0)
        std_moment1_x = np.sqrt(var_moment1_x)
        var_moment1_y = np.sum(np.square(np.multiply(weighted_var, np.subtract(y, moment1_y)))) / np.square(moment0)
        std_moment1_y = np.sqrt(var_moment1_y)
    except TypeError:
        std_moment1_x = None
        std_moment1_y = None

    return moment0, std_moment0, moment1_x, std_moment1_x, moment1_y, std_moment1_y


def gaussian_weights(shape, fwhm, center=None, normalization=1):
    if center is None:
        center = (shape[0]-1) / 2, (shape[1]-1) / 2

    y, x = np.mgrid[:shape[0], :shape[1]]
    r_squared = (y - center[0])**2 + (x - center[1])**2

    var = fwhm / 2.35
    weights = np.exp(-np.divide(r_squared, 2 * var)) / np.sqrt(2 * np.pi * var)

    return normalize_weights(weights=weights, normalization=normalization)


def uniform_weights(shape, normalization=1):
    return normalize_weights(np.ones(shape=shape), normalization=normalization)


def normalize_weights(weights, normalization=None):
    try:
        return weights * np.divide(normalization, np.sum(weights))
    except TypeError:
        return weights
