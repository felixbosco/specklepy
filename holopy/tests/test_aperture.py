import unittest
import numpy as np
from holopy.core.aperture import Aperture
from holopy.utils.imshow import imshow


class TestAperture(unittest.TestCase):

    def setUp(self):
        self.test_data = np.arange(0, 16*16).reshape((16, 16))

    def test_init(self):
        Aperture(8, 8, radius=4, data=self.test_data)

    def test_call(self):
        test_aperture = Aperture(8, 8, radius=4, data=self.test_data)
        test_aperture()
        # imshow(test_aperture())

if __name__ == "__main__":
    unittest.main()
