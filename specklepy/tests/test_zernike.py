import unittest
import numpy as np

from holopy.utils.zernike import Zernike
from holopy.utils.plot import imshow


class TestZernike(unittest.TestCase):

    def setUp(self):
        self.size = 256

    def test_init(self):
        z = Zernike()
        imshow(z.init_rho(self.size), title='Radius')
        imshow(z.init_phi(self.size), title='Azimuth')

    def test_init_from_vector(self):
        z = Zernike()
        coeffs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        # coeffs = np.random.rand((10))
        out = z(coeffs, size=128)
        imshow(out, title='Zernike polynomial {}'.format(coeffs))

    def test_init_from_keyword(self):
        z = Zernike()
        imshow(z.defocus(-1, 256), title='defocus')

if __name__ == "__main__":
    unittest.main()
