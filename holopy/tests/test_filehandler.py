import unittest
import numpy as np
from holopy.io.filehandler import FileHandler
from holopy.io.outfile import Outfile


class TestAperture(unittest.TestCase):

    def setUp(self):
        self.filename = "data/test/example_cube.fits"

    def test_init(self):
        FileHandler(self.filename)

    # def test_create_outfile(self):
    #     fh = FileHandler(self.filename)
    #     fh.create_outfile()

if __name__ == "__main__":
    unittest.main()
