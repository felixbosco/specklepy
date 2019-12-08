import unittest
import numpy as np

from specklepy.core.combine import weighted_combine
from specklepy.utils.plot import imshow



class TestCombine(unittest.TestCase):

    def setUp(self):
        self.shape = (4, 32, 32)
        self.data = np.random.rand(*self.shape)
        self.data[0] *= 3
        self.data[1] *= 1
        self.data[2] *= 4
        self.data[3] *= 5
        self.vars = np.random.rand(*self.shape) * 0.5


    def test_weighted_combin(self):
        mean, weights = weighted_combine(self.data, axis=0)
        mean, weights = weighted_combine(self.data, axis=0, vars=self.vars)
        imshow(mean)
        imshow(weights)



if __name__ == "__main__":
    unittest.main()
