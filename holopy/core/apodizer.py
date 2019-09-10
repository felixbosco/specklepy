import numpy as np
from astropy.modeling.models import Gaussian2D
from astropy.modeling.models import AiryDisk2D

class Apodizer(object):

    def __init__(self, function, **kwargs):

        if function == 'Gaussian':
            self.model = Gaussian2D(**kwargs)
        elif function == 'Airy':
            self.model = AiryDisk2D(**kwargs)
        else:
            raise ValueError("Function value <{}> for Apodizer class is not recognized!".format(function))

    def __call__(self, size, **kargs):
        y, x = np.mgrid[0:size, 0:size]
        return self.model(x, y, **kwargs)
