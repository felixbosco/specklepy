import numpy as np

from holopy.logging import logging

class Aperture(object):

    def __init__(self, x0, y0, width, data, mask='circular'):
        self.x0 = x0
        self.y0 = y0
        self.width = width

    @property
    def radius(self):
        return self.width
