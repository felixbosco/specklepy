import numpy as np
from astropy.modeling.models import Gaussian2D
from astropy.modeling.models import AiryDisk2D

from specklepy.exceptions import SpecklepyValueError
from specklepy.logging import logger
from specklepy.utils.transferfunctions import otf


def apodize(fourier_object, function='gaussian', crop=False, **kwargs):
    """Apodize an object with a kernel.

    Args:
        fourier_object (np.ndarray):
            Object to apodize.
        function (str, optional):
            Model function to initialize the apodization kernel along. Models are astropy.modelling.models Gaussian2D or
            AiryDisk2D. Default is 'gaussian'.
        crop (bool, optional):
            Crop corners of the PSF and set them to zero.
        kwargs (optional):
            kwargs are passed to the model function.
    """

    # Assert object type and shape
    if not isinstance(fourier_object, np.ndarray) and fourier_object.ndim != 2:
        raise ValueError("object must be a 2-dimensional np.ndarray!")
    if fourier_object.shape[0] != fourier_object.shape[1]:
        logger.warn("specklepy.core.apodization.apodize received a non quadratic input image. This may cause some "
                    "unpredictable results!")

    # Interpret function input and initialize model function
    center = ((fourier_object.shape[0] + 1) / 2, (fourier_object.shape[1] + 1) / 2)

    if 'gauss' in function.lower():
        if 'radius' in kwargs:
            # Copying the radius keyword argument into the proper Gaussian2D keywords
            kwargs['x_stddev'] = kwargs['radius']
            kwargs['y_stddev'] = kwargs['radius']
            del kwargs['radius']
        model = Gaussian2D(x_mean=center[0], y_mean=center[1], **kwargs)

    elif 'airy' in function.lower():
        model = AiryDisk2D(x_0=center[0], y_0=center[1], **kwargs)

    else:
        raise SpecklepyValueError("func apodize", "function", function, "either gaussian or airy")

    # Evaluate apodization function on object image grid and multiply on object
    y, x = np.mgrid[0:fourier_object.shape[0], 0:fourier_object.shape[1]]
    apodization_psf = model(x, y)

    # Crop corners of the PSF
    if crop:
        threshold = apodization_psf[0, int(center[1])]
        apodization_psf -= threshold
        apodization_psf = np.maximum(apodization_psf, 0.0)

    # Normalize to unity
    apodization_psf /= np.sum(apodization_psf)

    # Transform into Fourier space
    apodization_otf = otf(apodization_psf)
    return np.multiply(fourier_object, apodization_otf)
