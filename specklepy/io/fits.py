import os
import sys

from astropy.io import fits


def get_data(file_name, path=None, extension=None):

    # Construct path
    if path is None:
        path = file_name
    else:
        path = os.path.join(path, file_name)

    # Load data
    try:
        data = fits.getdata(path, extension)
    except (FileNotFoundError, KeyError) as e:
        sys.tracebacklimit = 0
        raise e

    # Return
    return data
