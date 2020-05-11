import unittest
import numpy as np

from specklepy.core.segmentation import Segmentation
from specklepy.utils.plot import imshow


class TestSegmentation(unittest.TestCase):

    def setUp(self):
        self.shape = (4, 32, 32)
        self.data = np.random.rand(*self.shape)
        self.data[0] *= 3
        self.data[1] *= 1
        self.data[2] *= 4
        self.data[3] *= 5
        self.vars = np.random.rand(*self.shape) * 0.5

    def test_init(self):
        seg = Segmentation(2, 2, self.shape)


if __name__ == "__main__":
    unittest.main()