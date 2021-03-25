import numpy as np


def first_moment_2d(array, var=None, weights=None):

    # Digest weights
    if weights is None:
        weights = np.ones(array.shape)
    else:
        raise NotImplementedError("Using weights for computing the first moment is not implemented yet!")
    weighted_array = np.multiply(weights, array)

    # Digest var
    if var is None:
        var = np.zeros(array.shape)
    weighted_var = np.multiply(weights, var)

    # Initialize arrays
    x, y = np.mgrid[:array.shape[0], :array.shape[1]]

    # Compute zeroth moment
    moment0 = np.sum(weighted_array)
    var_moment0 = np.sum(np.square(weighted_var))

    # Compute first moments
    moment1_x = np.average(x, weights=weighted_array, axis=(0, 1))
    moment1_y = np.average(y, weights=weighted_array, axis=(0, 1))

    # Compute uncertainties
    var_moment1_x = np.sum(np.square(np.multiply(weighted_var, np.subtract(x, moment1_x)))) / np.square(moment0)
    var_moment1_y = np.sum(np.square(np.multiply(weighted_var, np.subtract(y, moment1_y)))) / np.square(moment0)

    return moment0, np.sqrt(var_moment0), moment1_x, moment1_y, np.sqrt(var_moment1_x), np.sqrt(var_moment1_y)
