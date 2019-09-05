import unittest
from holopy.io.filehandler import FileHandler


class TestFileHandler(unittest.TestCase):

    def setUp(self):
        self.filename = "data/test/example_cube.fits"

    def test_init(self):
        FileHandler(self.filename)

    def test_str(self):
        fh = FileHandler(self.filename)
        print(fh)

if __name__ == "__main__":
    unittest.main()
