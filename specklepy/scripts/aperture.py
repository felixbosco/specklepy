import numpy as np
import os

from astropy.io import fits
from astropy.table import Table, Column

from specklepy import imshow
from specklepy.core import Aperture
from specklepy.logging import logger


def aperture_analysis(file, index, radius, out_file=None, pixel_scale=1, recenter=False, debug=False):

    if out_file is None:
        out_file = 'aperture_' + os.path.basename(file).replace(".fits", ".dat")

    # Initialize the aperture
    aperture = Aperture(index, radius, file_name=file, crop=not recenter)

    # Recenter aperture on peak
    if recenter:
        aperture.center_on_peak()
        aperture.crop()

    # Initialize the output table
    out_table = Table()

    # PSF profile analysis
    radius, flux, flux_err = aperture.get_psf_profile()
    out_table.add_column(Column(data=radius, name='Radius'))
    out_table.add_column(Column(data=flux, name='Flux'))
    out_table.add_column(Column(data=flux_err, name='FluxError'))

    # Power spectrum analysis
    try:
        radius, mean, std = aperture.get_power_spectrum_profile()
        spat_freq = aperture.spatial_frequency(pixel_scale=pixel_scale)
        spat_wave = aperture.spatial_wavelength(pixel_scale=pixel_scale)
        out_table.add_column(Column(data=spat_freq, name='SpatialFrequency'))
        out_table.add_column(Column(data=spat_wave, name='SpatialWavelength'))
        out_table.add_column(Column(data=mean, name='AmplitudeMean'))
        out_table.add_column(Column(data=std, name='AmplitudeStd'))
    except IndexError:
        logger.error("Image data is not a cube. Skipping time evolution")

    # Store output table
    logger.info(f"Storing table to file {out_file!r}")
    out_table.write(out_file, overwrite=True, format='ascii.fixed_width')


def get_psf_1d(file, index, radius, out_file=None, normalize=None, debug=False):

    if isinstance(index, list):
        if len(index) == 1:
            if index[0] == 0:
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
    aperture = Aperture(index, radius, file_name=file, crop=True)
    if debug:
        imshow(aperture.get_integrated(), maximize=False)

    # Extract PSF profile
    logger.info(f"Extracting PSF profile from file {file}")
    x_data, y_data, e_data = aperture.get_psf_profile()

    # Normalize profile
    if normalize == 'peak':
        y_data /= y_data[0]
        e_data /= y_data[0]
    elif normalize == 'aperture':
        y_data /= y_data[-1]
        e_data /= y_data[-1]
    elif normalize is not None:
        raise ValueError("Normalize must be either 'peak', 'aperture, or None!'")

    # Save encircled energy data to outfile
    out_table = Table(data=[x_data, y_data, e_data], names=['Radius', 'Flux', 'dFlux'])
    logger.info(f"Store PSF profile to {out_file}")
    out_table.write(out_file, overwrite=True, format='ascii.fixed_width')


def get_psf_variation(file, index, radius, out_file=None, normalize=None, debug=False):
    if isinstance(index, list):
        if len(index) == 1:
            if index[0] == 0:
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
    aperture = Aperture(index, radius, file_name=file, crop=True)
    if debug:
        imshow(aperture.get_integrated(), maximize=False)

    # Extract PSF profile
    logger.info(f"Extracting PSF profile from file {file}")
    x_data, y_data, e_data = aperture.get_psf_variance()

    # Normalize profile
    if normalize == 'peak':
        y_data /= y_data[0]
        e_data /= y_data[0]
    elif normalize == 'aperture':
        y_data /= y_data[-1]
        e_data /= y_data[-1]
    elif normalize is not None:
        raise ValueError("Normalize must be either 'peak', 'aperture, or None!'")

    # Save encircled energy data to outfile
    out_table = Table(data=[x_data, y_data, e_data], names=['Radius', 'Variance', 'dVariance'])
    logger.info(f"Store PSF profile to {out_file}")
    out_table.write(out_file, overwrite=True, format='ascii.fixed_width')
