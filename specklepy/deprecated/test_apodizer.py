import unittest
from matplotlib.colors import LogNorm
from holopy.core.apodizer import Apodizer
from holopy.utils.plot import imshow


class TestApodizer(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        Apodizer('Gaussian', 256)
        Apodizer('Airy', 256)
        with self.assertRaises(ValueError):
            Apodizer('Nonsense', 256)

    def test_call(self):
        apodizer = Apodizer('Gaussian', 256, radius=10)
        imshow(apodizer())

if __name__ == "__main__":
    unittest.main()
