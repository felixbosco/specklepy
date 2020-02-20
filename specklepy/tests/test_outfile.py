import unittest
import numpy as np
from specklepy.io.filemanager import FileManager
from specklepy.io.outfile import Outfile


class TestOutfile(unittest.TestCase):

    def setUp(self):
        self.path = "specklepy/tests/files/"
        self.file = "test_outfile.fits"

    def test_init(self):
        with self.assertRaises(RuntimeError):
            Outfile(None)
        Outfile(filename=self.path+self.file, shape=None, cards={"RECONSTRUCTION": "Test"})
        Outfile(filename=self.path+self.file, shape=(10, 10), extensions={'name': 'var', 'shape': (10, 20)}, cards={"RECONSTRUCTION": "Test"}, timestamp=True)

    def test_set_data(self):
        pass

    def test_new_extension(self):
        outfile = Outfile(filename=self.path+self.file, shape=(10, 10), cards={"RECONSTRUCTION": "Test"})
        outfile.new_extension(name='VAR', data=np.zeros((10,10)))


if __name__ == "__main__":
    unittest.main()
