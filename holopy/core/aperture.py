import numpy as np
from copy import copy

# from holopy.logging import logging

class Aperture(object):

    def __init__(self, x0, y0, radius, data, mask='circular', subset_only=True):
        self.x0 = x0
        self.y0 = y0
        self.radius = radius

        if subset_only:
            self.data = copy(data[x0 - radius : x0 + radius + 1, y0 - radius : y0 + radius + 1])
        else:
            self.data = copy(data)

        if mask is None:
            pass
        elif mask is 'circular':
            self.data = np.ma.masked_array(self.data, mask=self.make_mask())
        else:
            raise ValueError("Mask type '{}' of Aperture instance not understood.".format(mask))

    @property
    def width(self):
        return self.radius * 2 + 1

    def __call__(self):
        return self.data

    def make_mask(self):
        xx, yy = np.mgrid[:self.data.shape[0], :self.data.shape[1]]
        distance_map = np.sqrt(np.square(xx - self.radius) + np.square(yy - self.radius))
        return np.ma.masked_greater(distance_map, self.radius).mask
