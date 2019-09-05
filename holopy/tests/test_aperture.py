import unittest
import numpy as np
from holopy.core.aperture import Aperture
from holopy.utils.imshow import imshow


class TestAperture(unittest.TestCase):

    def setUp(self):
        size = 16
        self.test_data = np.arange(0, size*size).reshape((size, size))
        size_large = 128
        self.test_data_large = np.arange(0, size_large*size_large).reshape((size_large, size_large))


    def test_init(self):
        Aperture(8, 8, radius=4, data=self.test_data)

    def test_call(self):
        test_aperture = Aperture(8, 8, radius=4, data=self.test_data)
        test_aperture()
        test_aperture = Aperture(64, 64, radius=16, data=self.test_data_large)
        test_aperture()
        imshow(test_aperture())

if __name__ == "__main__":
    unittest.main()