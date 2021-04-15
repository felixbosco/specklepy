import unittest
import numpy as np

from specklepy.core.aperture import Aperture
from specklepy.plotting import utils


class TestTransferFunctions(unittest.TestCase):

    def setUp(self):
        size = 256
        amp = np.ones((size, size))
        amp = Aperture(128, 128, 64, file_name=amp, crop=False)().filled(0)
        pha = np.random.rand(size, size) * 2*np.pi * 2**(-1)
        self.aperture = amp * np.exp(1j * pha)

    def test_desaturate_color(self):
        utils.desaturate_color('tab:blue', number_colors=3)
        utils.desaturate_color('#00FF00', number_colors=3)
        utils.desaturate_color((0, 1, 0), number_colors=3)


if __name__ == "__main__":
    unittest.main()
