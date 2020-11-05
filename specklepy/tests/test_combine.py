import unittest
import numpy as np

from specklepy.utils.combine import weighted_mean, time_difference
from specklepy.plotting.plot import imshow


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
        mean, weights = weighted_mean(self.data, axis=0)
        mean, weights = weighted_mean(self.data, axis=0, vars=self.vars)
        imshow(mean)
        imshow(weights)

    def test_time_difference(self):
        t0 = '2020-04-23 13:02:31.397899'
        self.assertTrue(time_difference(t0, t0) == 0)
        times = ['2018-11-14 16:51:37.070953', '2020-04-17 21:49:47.877587',
                 '2020-04-23 13:02:31.397899', '2020-04-23 13:14:11.297996']
        deltas = time_difference(t0, times)
        self.assertTrue(type(deltas) is np.ndarray)


if __name__ == "__main__":
    unittest.main()
