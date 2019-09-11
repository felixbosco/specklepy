import unittest
import numpy as np
from matplotlib.colors import LogNorm

from holopy.core.aperture import Aperture
from holopy.utils.transferfunctions import psf, otf, mtf, powerspec
from holopy.utils.imshow import imshow


class TestTransferFunctions(unittest.TestCase):

    def setUp(self):
        size = 256
        data = np.ones((size, size))
        self.aperture = Aperture(x0=128, y0=128, radius=32, data=data, subset_only=False)().filled(0)
        imshow(self.aperture, title="Aperture")

    def test_psf(self):
        psf = psf(self.aperture)
        imshow(psf, title="PSF", norm=LogNorm())

if __name__ == "__main__":
    unittest.main()
