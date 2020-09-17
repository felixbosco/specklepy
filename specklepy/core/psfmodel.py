import numpy as np

from astropy.modeling.models import Gaussian2D, AiryDisk2D

from specklepy.exceptions import SpecklepyValueError


class PSFModel(object):

    def __init__(self, type, shape=None, **kwargs):

        # Store input
        self.type = type
        self.shape = shape

        # Derive center of model shape
        if self.shape is None:
            self.center = (0, 0)
        else:
            self.center = ((self.shape[0] + 1) / 2, (self.shape[1] + 1) / 2)

        # Initialize model
        if 'gauss' in self.type.lower():
            if 'radius' in kwargs:
                # Copying the radius keyword argument into the proper Gaussian2D keywords
                kwargs['x_stddev'] = kwargs['radius']
                kwargs['y_stddev'] = kwargs['radius']
                del kwargs['radius']
            self.model = Gaussian2D(x_mean=self.center[0], y_mean=self.center[1], **kwargs)

        elif 'airy' in self.type.lower():
            self.model = AiryDisk2D(x_0=self.center[0], y_0=self.center[1], **kwargs)

        else:
            raise SpecklepyValueError('PSFModel', argname='type', argvalue=type, expected="either 'Gaussian' or 'Airy'")

    def __call__(self, shape):

        # Update shape and center according to input shape
        self.shape = shape
        self.center = ((self.shape[0] + 1) / 2, (self.shape[1] + 1) / 2)

        # Update self.model center coordinates
        from IPython import embed
        if isinstance(self.model, Gaussian2D):
            self.model.x_mean, self.model.y_mean = list(self.center)
        elif isinstance(self.model, AiryDisk2D):
            self.model.x_0, self.model.y_0 = list(self.center)

        # Create coordinate grid
        y, x = np.mgrid[0:self.shape[0], 0:self.shape[1]]

        # Evaluate model on grid and return
        return self.model(x, y)
