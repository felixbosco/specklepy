import unittest
import numpy as np
from specklepy.io.filemanager import FileManager
from specklepy.io.outfile import Outfile


class TestOutfile(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        with self.assertRaises(RuntimeError):
            Outfile(None)
        Outfile(filename="specklepy/tests/files/test_outfile.fits", shape=None, cards={"RECONSTRUCTION": "Test"})
        Outfile(filename="specklepy/tests/files/test_outfile.fits", shape=(10, 10), extensions='var', cards={"RECONSTRUCTION": "Test"}, timestamp=True)

    def test_set_data(self):
        pass

if __name__ == "__main__":
    unittest.main()
