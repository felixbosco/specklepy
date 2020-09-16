import numpy as np
import os

from astropy.table import Table

from specklepy.core.aperture import Aperture
from specklepy.utils.plot import imshow


def get_psf_1d(file, index, radius, out_file=None, normalize=None, debug=False):
    index = tuple(index)

    if file is None:
        raise RuntimeError("No file was provided!")

    if out_file is None:
        out_file = "psf_" + os.path.basename(file).replace(".fits", ".dat")

    # Initialize the aperture
    aperture = Aperture(index, radius, data=file, crop=True)
    if debug:
        imshow(aperture.get_integrated(), maximize=False)

    xdata, ydata = aperture.get_psf_profile()

    if normalize == 'peak':
        ydata /= ydata[0]
    elif normalize == 'aperture':
        ydata /= ydata[-1]
    elif normalize is not None:
        raise ValueError("Normalize must be either 'peak', 'aperture, or None!'")

    # Save encircled energy data to outfile
    out_table = Table(data=[xdata, ydata], names=['Radius', 'Flux'])
    out_table.write(out_file, overwrite=True, format='ascii.fixed_width')
