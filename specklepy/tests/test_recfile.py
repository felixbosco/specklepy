import unittest
import numpy as np
from specklepy.io.filemanager import FileManager
from specklepy.io.recfile import RECfile


class TestRECfile(unittest.TestCase):

    def setUp(self):
        self.FileManager = FileManager("data/test/example_cube.fits")

    def test_init(self):
        RECfile(filename="data/test/test_recfile.fits", files=self.FileManager.files, cards={"RECONSTRUCTION": "Test"})

    def test_set_data(self):
        pass

if __name__ == "__main__":
    unittest.main()
