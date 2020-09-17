from specklepy.logging import logger
from specklepy.utils.transforms import rad2mas


def first_airy_zero(wavelength, diameter):
    return 1.22 * wavelength / diameter


def first_airy_zero_to_gaussian_sigma(val):
    return val / 1.22 * 0.42


def get_resolution_parameters(wavelength, diameter, pixel_scale=None):
    airy_zero = first_airy_zero(wavelength=wavelength, diameter=diameter)
    logger.info(f"First airy zero: {airy_zero:.3e} rad | {rad2mas(airy_zero):.3f} mas")
    if pixel_scale:
        logger.info(f"First airy zero: {rad2mas(airy_zero) / pixel_scale:.3f} pix")

    sigma = airy_zero / 1.22 * 0.42
    logger.info(f"Gaussian sigma: {sigma:.3e} rad | {rad2mas(sigma):.3f} mas")
    if pixel_scale:
        logger.info(f"Gaussian sigma: {rad2mas(sigma) / pixel_scale:.3f} pix")
