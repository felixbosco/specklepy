import unittest
import numpy as np
from specklepy.core.aperture import Aperture
from specklepy.plotting.utils import imshow


class TestAperture(unittest.TestCase):

    def setUp(self):
        size = 16
        self.test_data = np.arange(0, size*size).reshape((size, size))
        size_large = 128
        self.test_data_large = np.arange(0, size_large*size_large).reshape((size_large, size_large))

    def test_init(self):
        Aperture(8, 8, 4, file_name=self.test_data)
        # Non-integer indizes:
        Aperture(8.2, 8.3, 4, file_name=self.test_data)
        # Non-interger radius:
        with self.assertRaises(ValueError):
            Aperture(8.2, 8.3, 3.7, file_name=self.test_data)

    def test_crop(self):
        aperture = Aperture(8, 8, 4, file_name=self.test_data, crop=True)
        aperture.crop()
        # imshow(aperture())
        aperture = Aperture(8, 8, 4, file_name=self.test_data, crop=False)
        # imshow(aperture())
        aperture.crop()
        # imshow(aperture())

    def test_encircled_energy(self):
        aperture = Aperture(8, 8, 4, file_name=self.test_data)
        aperture.get_encircled_energy(save_to='specklepy/tests/files/analysis/example_encircled_energy.dat')

    def test_call(self):
        imshow(self.test_data)
        test_aperture = Aperture(8, 8, 4, file_name=self.test_data)
        imshow(test_aperture.data)
        test_aperture = Aperture(64, 64, 16, file_name=self.test_data_large)
        imshow(test_aperture.data)


if __name__ == "__main__":
    unittest.main()
