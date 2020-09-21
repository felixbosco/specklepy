import numpy as np
import os

from astropy.io import fits
from astropy.table import Table

from specklepy.core.aperture import Aperture
from specklepy.logging import logger
from specklepy.utils.plot import imshow


def get_psf_1d(file, index, radius, out_file=None, normalize=None, debug=False):

    if isinstance(index, list):
        if len(index) is 1:
            if index[0] is 0:
                logger.info(f"Estimate image intensity peak and use as aperture index")
                image = fits.getdata(file)
                index = np.unravel_index(np.argmax(image), image.shape)
                logger.info(f"Index is set to {index}")
            else:
                index = (index[0], index[0])
        index = tuple(index)

    if file is None:
        raise RuntimeError("No file was provided!")

    if out_file is None:
        out_file = "psf_" + os.path.basename(file).replace(".fits", ".dat")

    # Initialize the aperture
    aperture = Aperture(index, radius, data=file, crop=True)
    if debug:
        imshow(aperture.get_integrated(), maximize=False)

    # Extract PSF profile
    logger.info(f"Extracting PSF profile from file {file}")
    xdata, ydata, edata = aperture.get_psf_profile()

    # Normalize profile
    if normalize == 'peak':
        ydata /= ydata[0]
        edata /= ydata[0]
    elif normalize == 'aperture':
        ydata /= ydata[-1]
        edata /= ydata[-1]
    elif normalize is not None:
        raise ValueError("Normalize must be either 'peak', 'aperture, or None!'")

    # Save encircled energy data to outfile
    out_table = Table(data=[xdata, ydata, edata], names=['Radius', 'Flux', 'dFlux'])
    logger.info(f"Store PSF profile to {out_file}")
    out_table.write(out_file, overwrite=True, format='ascii.fixed_width')


def get_psf_variation(file, index, radius, out_file=None, normalize=None, debug=False):
    if isinstance(index, list):
        if len(index) is 1:
            if index[0] is 0:
                logger.info(f"Estimate image intensity peak and use as aperture index")
                image = fits.getdata(file)
                if image.ndim == 3:
                    image = np.sum(image, axis=0)
                index = np.unravel_index(np.argmax(image), image.shape)
                logger.info(f"Index is set to {index}")
            else:
                index = (index[0], index[0])
        index = tuple(index)

    if file is None:
        raise RuntimeError("No file was provided!")

    if out_file is None:
        out_file = "var_" + os.path.basename(file).replace(".fits", ".dat")

    # Initialize the aperture
    aperture = Aperture(index, radius, data=file, crop=True)
    if debug:
        imshow(aperture.get_integrated(), maximize=False)

    # Extract PSF profile
    logger.info(f"Extracting PSF profile from file {file}")
    xdata, ydata, edata = aperture.get_psf_variance()

    # Normalize profile
    if normalize == 'peak':
        ydata /= ydata[0]
        edata /= ydata[0]
    elif normalize == 'aperture':
        ydata /= ydata[-1]
        edata /= ydata[-1]
    elif normalize is not None:
        raise ValueError("Normalize must be either 'peak', 'aperture, or None!'")

    # Save encircled energy data to outfile
    out_table = Table(data=[xdata, ydata, edata], names=['Radius', 'Variance', 'dVariance'])
    logger.info(f"Store PSF profile to {out_file}")
    out_table.write(out_file, overwrite=True, format='ascii.fixed_width')
