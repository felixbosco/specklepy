import unittest
import numpy as np
from matplotlib.colors import LogNorm

from specklepy.core.aperture import Aperture
from specklepy.utils.transferfunctions import psf, otf, mtf, powerspec
from specklepy.utils.plot import imshow, plot_powerspec1d


class TestTransferFunctions(unittest.TestCase):

    def setUp(self):
        size = 256
        amp = np.ones((size, size))
        amp = Aperture(x0=128, y0=128, radius=64, data=amp, subset_only=False)().filled(0)
        pha = np.random.rand(size, size) * 2*np.pi * 2**(-1)
        self.aperture = amp * np.exp(1j * pha)

    def test_plot_powerspec1d(self):
        psf_image = psf(self.aperture)
        imshow(psf_image, title="PSF", norm=LogNorm())
        plot_powerspec1d(psf_image)

if __name__ == "__main__":
    unittest.main()
