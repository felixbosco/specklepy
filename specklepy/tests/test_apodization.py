import unittest
import numpy as np
from specklepy.core.apodization import apodize
from specklepy.plotting.utils import imshow


class TestApodization(unittest.TestCase):

    def setUp(self):
        self.object = np.zeros((256, 256))
        self.object[128, 128] = 10

    def test_call(self):
        apodize(self.object, 'Gaussian', radius=16)
        apodize(self.object, 'Airy', radius=16)
        imshow(np.abs(apodize(self.object, 'Gaussian', radius=20)))

    def test_errors(self):


        with self.assertRaises(ValueError):
            apodize(self.object, 'Nonsense')

if __name__ == "__main__":
    unittest.main()
