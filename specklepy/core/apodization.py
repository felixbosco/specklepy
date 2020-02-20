import numpy as np
from astropy.modeling.models import Gaussian2D
from astropy.modeling.models import AiryDisk2D

from specklepy.utils.transferfunctions import otf
from specklepy.utils.plot import imshow
from specklepy.logging import logger



def apodize(object, function='gaussian', crop=False, **kwargs):
    """Apodize an object with a given model function.

    Long description...

    Args:
        object (np.ndarray): Object to apodize.
        function (str, optional): Function to initialize the apodization model
            along. Models are astropy.modelling.models Gaussian2D or AiryDisk2D.
            Default is 'gaussian'.
        kwargs (optional): kwargs are passed to the model function.
    """


    # Assert object type and shape
    if not isinstance(object, np.ndarray) and object.ndim != 2:
        raise ValueError("object must be a 2-dimensional np.ndarray!")
    if object.shape[0] != object.shape[1]:
        logger.warn("specklepy.core.apodization.apodize received a non quadratic input image. This may cause some unpredictable results!")



    # Interprete function input and initialize model function
    _function = function.lower()
    center = ((object.shape[0] + 1) / 2, (object.shape[1] + 1) / 2)

    if 'gauss' in _function:
        if 'radius' in kwargs:
            # Copying the radius keyword argument into the proper Gaussian2D keywords
            kwargs['x_stddev'] = kwargs['radius']
            kwargs['y_stddev'] = kwargs['radius']
            del kwargs['radius']
        model = Gaussian2D(x_mean=center[0], y_mean=center[1], **kwargs)

    elif 'airy' in _function:
        model = AiryDisk2D(x_0=center[0], y_0=center[1], **kwargs)

    else:
        raise ValueError("specklepy.core.apodization.apodize received unknown value for function '{}'!".format(function))

    # Evaluate apodization function on object image grid and multiply on object
    y, x = np.mgrid[0:object.shape[0], 0:object.shape[1]]
    apodization_psf = model(x, y)

    # Crop corners of the PSF
    if crop:
        threshold = apodization_psf[0, int(center[1])]
        apodization_psf -= threshold
        apodization_psf = np.maximum(apodization_psf, 0.0)

    # Transform into Fourier space
    apodization_otf = otf(apodization_psf)
    return np.multiply(object, apodization_otf)
