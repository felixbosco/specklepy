import numpy as np


def rad2deg(val):
    """Transform an angle from units of radians into degrees."""
    return np.rad2deg(val)


def rad2mas(val):
    """Transform an angle from units of radians into milliarcseconds."""
    return np.rad2deg(val) * 60**2 * 1000


def ang2pix(val, pixel_scale):
    """Transform an angle from angular units into pixels."""
    return val / pixel_scale


def pix2ang(val, pixel_scale):
    """Transform an angle from units of pixels into angular units."""
    return val * pixel_scale


def flux2mag(flux, flux_reference):
    """Transform a brightness value from units of flux into magnitudes."""
    return -2.5 * np.log10(np.divide(flux, flux_reference))


def mag2flux(mag, flux_reference=1):
    """Transform a brightness value from units of magnitudes into a flux, in the units of the reference flux."""
    return np.power(10, np.divide(mag, -2.5)) * flux_reference
