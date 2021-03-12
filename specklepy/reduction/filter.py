import numpy as np


def hot_pixel_mask(image, threshold=5):
    """Identify hot pixels via a 2nd derivative method and create a hot pixel mask with hot pixels `True`.

    Arguments:
        image (np.ndarray):
            Image for which the mask shall be created.
        threshold (int or float, optional):
            Sigma-clipping threshold, used to evaluate whether the 2nd derivative is too large.

    Returns:
        mask (np.ndarray):
            Boolean array of the shape of the image input. Hot pixels are `True`.
    """

    # Compute second derivatives
    ddx = np.diff(image, n=2)
    ddy = np.diff(image, n=2, axis=-2)

    # Pad derivatives to return to image size
    ddx = np.pad(ddx, ((0, 0), (1, 1)))
    ddy = np.pad(ddy, ((1, 1), (0, 0)))

    # Derive masks
    mask_x = np.abs(ddx) > threshold * np.std(ddx)
    mask_y = np.abs(ddy) > threshold * np.std(ddy)

    return mask_x | mask_y


def mask_hot_pixels(image, threshold=5):
    return np.ma.masked_array(image, mask=hot_pixel_mask(image=image, threshold=threshold))


def fill_hot_pixels(image, fill_value=0, threshold=5):
    return mask_hot_pixels(image=image, threshold=threshold).filled(fill_value)


def variable_pixel_mask(cube, var, threshold=5):
    std = np.std(cube, axis=0)
    return np.divide(std, np.sqrt(var)) > threshold


def mask_variable_pixels(cube, var, threshold=5):
    integrated = np.sum(cube, axis=0)
    return np.ma.masked_array(integrated, mask=variable_pixel_mask(cube=cube, var=var, threshold=threshold))


def fill_variable_pixels(cube, var, fill_value=0, threshold=5):
    return mask_variable_pixels(cube=cube, var=var, threshold=threshold).filled(fill_value)


def bad_pixel_mask(cube, var, threshold=5):
    hpm = hot_pixel_mask(np.sum(cube, axis=0), threshold=threshold)
    vpm = variable_pixel_mask(cube=cube, var=var, threshold=threshold)
    return np.logical_or(hpm, vpm)
