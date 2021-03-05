import unittest
import numpy as np
from matplotlib.colors import LogNorm

from specklepy.core.aperture import Aperture
from specklepy.utils.transferfunctions import psf
from specklepy.plotting.utils import imshow, plot_powerspec1d


class TestTransferFunctions(unittest.TestCase):

    def setUp(self):
        size = 256
        amp = np.ones((size, size))
        pha = np.random.rand(size, size) * 1e1
        # imshow(pha, title="Phase")
        amp = Aperture(128, 128, 64, data=amp, crop=False).data
        self.aperture = amp * np.exp(1j * pha)
        # imshow(np.abs(self.aperture), title="abs Aperture")

    def test_psf(self):
        psf_image = psf(self.aperture)
        imshow(psf_image, title="PSF", norm=LogNorm())

    def test_powerspec1d(self):
        psf_image = psf(self.aperture)
        plot_powerspec1d(psf_image, title='simple')
        plot_powerspec1d(psf_image, title='no average', average=False)
        # plot_powerspec1d(psf_image, title='pixel_scale', pixel_scale=15*u.mas)

if __name__ == "__main__":
    unittest.main()
