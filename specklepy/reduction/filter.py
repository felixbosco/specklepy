import numpy as np


def mask_hot_pixels(image, sigma_threshold=3):
    # Compute second derivatives
    ddx = np.diff(image, n=2)
    ddy = np.diff(image, n=2, axis=-2)

    # Pad derivatives to return to image size
    ddx = np.pad(ddx, ((0, 0), (1, 1)))
    ddy = np.pad(ddy, ((1, 1), (0, 0)))

    # Derive masks
    mask_x = np.abs(ddx) > sigma_threshold * np.std(ddx)
    mask_y = np.abs(ddy) > sigma_threshold * np.std(ddy)

    return np.ma.masked_array(image, mask=mask_x | mask_y)


def hot_pixel_mask(image, sigma_threshold=5):
    return mask_hot_pixels(image=image, sigma_threshold=sigma_threshold).mask
