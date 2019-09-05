import unittest
import numpy as np
from holopy.io.filehandler import FileHandler
from holopy.io.outfile import Outfile


class TestOutfile(unittest.TestCase):

    def setUp(self):
        self.FileHandler = FileHandler("data/test/example_cube.fits")

    def test_init(self):
        Outfile(file_list=self.FileHandler.files, filename="data/test/test_outfile.fits")

    def test_set_data(self):
        pass

if __name__ == "__main__":
    unittest.main()
