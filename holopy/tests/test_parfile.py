import unittest
from holopy.io.parfile import Parfile


class TestParfile(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        Parfile(filename="data/test/test_parfile.txt")

if __name__ == "__main__":
    unittest.main()
