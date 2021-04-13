import numpy as np

from specklepy.core.psfmodel import PSFModel
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

    # Interpret function input and comput apodization PSF
    psf_model = PSFModel(type=function, **kwargs)
    apodization_psf = psf_model(fourier_object.shape)

    # Crop corners of the PSF
    if crop:
        threshold = apodization_psf[0, int(psf_model.center[1])]
        apodization_psf -= threshold
        apodization_psf = np.maximum(apodization_psf, 0.0)

    # Normalize to unity
    apodization_psf /= np.sum(apodization_psf)

    # Transform into Fourier space
    apodization_otf = otf(apodization_psf)
    return np.multiply(fourier_object, apodization_otf)
