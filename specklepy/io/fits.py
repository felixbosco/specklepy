import os
import sys

from astropy.io import fits


def get_data(file_name, path=None, extension=None, squeeze=True, dtype=None, ignore_missing_extension=False):

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

    # Squeeze array
    if squeeze:
        data = data.squeeze()

    # Cast data type if requested
    if dtype is not None:
        data = data.astype(dtype=dtype)

    # Return
    return data
