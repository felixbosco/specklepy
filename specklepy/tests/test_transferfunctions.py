import unittest
import numpy as np
from matplotlib.colors import LogNorm

from specklepy.core.aperture import Aperture
from specklepy.utils.transferfunctions import psf
from specklepy.plotting.utils import imshow


class TestTransferFunctions(unittest.TestCase):

    def setUp(self):
        size = 256
        amp = np.ones((size, size))
        pha = np.random.rand(size, size) * 1e1
        amp = Aperture(128, 128, 64, file_name=amp, crop=False).data
        self.aperture = amp * np.exp(1j * pha)

    def test_psf(self):
        psf_image = psf(self.aperture)
        imshow(psf_image, title="PSF", norm=LogNorm())


if __name__ == "__main__":
    unittest.main()
