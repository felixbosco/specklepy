import numpy as np
from astropy.modeling.models import Gaussian2D
from astropy.modeling.models import AiryDisk2D

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


    def apodize(self):
        y, x = np.mgrid[0:self.size, 0:self.size]
        OTF = (np.fft.fft2(self.model(x, y)))


    def PSF(self, size=None):
        if size is not None:
            self.size = size
        y, x = np.mgrid[0:self.size, 0:self.size]
        self.psf = self.model(x, y)


    def OTF
