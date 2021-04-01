import os
import sys

from astropy.io import fits


def get_data(file_name, path=None, extension=None, ignore_missing_extension=False):

    # Construct path
    if path is None:
        path = file_name
    else:
        path = os.path.join(path, file_name)

    # Load data
    try:
        data = fits.getdata(path, extension)
    except FileNotFoundError as e:
        sys.tracebacklimit = 0
        raise e
    except KeyError as e:
        if ignore_missing_extension:
            return None
        sys.tracebacklimit = 0
        raise e

    # Return
    return data
