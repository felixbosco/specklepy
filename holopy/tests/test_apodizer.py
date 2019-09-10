import unittest
from holopy.core.apodizer import Apodizer


class TestApodizer(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        Apodizer('Gaussian')
        Apodizer('Airy')
        with self.assertRaises(ValueError):
            Apodizer('Nonsense')

    def test_call(self):
        pass

if __name__ == "__main__":
    unittest.main()
