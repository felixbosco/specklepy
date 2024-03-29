import os
import sys

from astropy.io import fits

from specklepy.utils import array


def get_data(file_name, path=None, extension=None, squeeze=True, dtype=None, ignore_missing_extension=False):
    """Read data from a FITS file.

    Args:
        file_name (str):
            Name of the file to read data from.
        path (str, optional):
            The path will be joined with `file_name`, if provided.
        extension (str, optional):
            Name of a FITS file extension. The behavior for missing extensions depends on the
            `ignore_missing_extension` argument.
        squeeze (bool, optional):
            If requested, the data array will be squeezed, i.e. length-1 dimensions are dropped.
        dtype (data type, optional):
            Data will be casted to this data type, if provided.
        ignore_missing_extension (bool, optional):
            Set the behavior for missing extensions. If `False`, a KeyError is raised, and otherwise the function
            returns a `None`.

    Returns:
        data (np.ndarray or None-type):
            Data array from the requested FITS file extension. Can return `None`, if extension is not contained.
    """

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


def get_header(file_name, path=None, extension=None, card=None, ignore_missing_extension=False):
    """Read data from a FITS file.

    Args:
        file_name (str):
            Name of the file to read data from.
        path (str, optional):
            The path will be joined with `file_name`, if provided.
        extension (str, optional):
            Name of a FITS file extension. The behavior for missing extensions depends on the
            `ignore_missing_extension` argument.
        card (str, optional):
            Name of a specific header card.
        ignore_missing_extension (bool, optional):
            Set the behavior for missing extensions. If `False`, a KeyError is raised, and otherwise the function
            returns a `None`.

    Returns:
        header (astropy.io.fits.Header or None-type):
            Header of the requested FITS file extension. Can return `None`, if extension is not contained in the file.
        value (any):
            Value of a specific card, if `card` is not None.
    """

    # Construct path
    if path is None:
        path = file_name
    else:
        path = os.path.join(path, file_name)

    # Load data
    try:
        header = fits.getheader(path, extension)
    except FileNotFoundError as e:
        sys.tracebacklimit = 0
        raise e
    except KeyError as e:
        if ignore_missing_extension:
            return None
        sys.tracebacklimit = 0
        raise e

    # Return value of requested card
    if card is not None:
        return header[card]

    # Return
    return header


def get_frame_number(file, path=None, extension=None):
    data = get_data(file_name=file, path=path, extension=extension)
    return array.frame_number(data)


def get_frame_shape(file, path=None, extension=None):
    data = get_data(file_name=file, path=path, extension=extension)
    return array.frame_shape(data)

