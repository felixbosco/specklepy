import unittest
import numpy as np
from holopy.io.filemanager import FileManager
from holopy.io.outfile import Outfile


class TestOutfile(unittest.TestCase):

    def setUp(self):
        self.FileManager = FileManager("data/test/example_cube.fits")

    def test_init(self):
        Outfile(files=self.FileManager.files, filename="data/test/test_outfile.fits", cards={"RECONSTRUCTION": "Test"})

    def test_set_data(self):
        pass

if __name__ == "__main__":
    unittest.main()
