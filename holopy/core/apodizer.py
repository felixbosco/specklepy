import numpy as np
from astropy.modeling.models import Gaussian2D
from astropy.modeling.models import AiryDisk2D

from holopy.utils.transferfunctions import otf

class Apodizer(object):

    def __init__(self, function, size, **kwargs):

        self.size = size
        self.center = (self.size + 1) / 2

        if function == 'Gaussian':
            if 'radius' in kwargs:
                # Copying the radius keyword argument into the proper Gaussian2D keywords
                kwargs['x_stddev'] = kwargs['radius']
                kwargs['y_stddev'] = kwargs['radius']
                del kwargs['radius']
            self.model = Gaussian2D(x_mean=self.center, y_mean=self.center, **kwargs)
        elif function == 'Airy':
            self.model = AiryDisk2D(x_0=self.center, y_0=self.center, **kwargs)
        else:
            raise ValueError("Function value <{}> for Apodizer class is not recognized!".format(function))


    def __call__(self):
        y, x = np.mgrid[0:self.size, 0:self.size]
        return self.model(x, y)


    def apodize(self, object):
        y, x = np.mgrid[0:self.size, 0:self.size]
        apodization_psf = self.model(x, y)
        apodization_otf = otf(apodization_psf)
        return np.multiply(object, apodization_otf)
        # tmp = np.multiply(object, apodization_otf)
        # return otf(tmp, inverse=True)
