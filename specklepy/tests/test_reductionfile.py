import unittest
import os
import numpy as np
from specklepy.io.filemanager import FileManager
from specklepy.io.reductionfile import ReductionFile


class TestReductionFile(unittest.TestCase):

    def setUp(self):
        self.path = "specklepy/tests/files/"
        self.parentFile = os.path.join(self.path, "reduction/MasterFlat.fits")

    def test_init(self):
        ReductionFile(self.parentFile, prefix='rtest_', path=os.path.join(self.path, 'reduction/'), reduction='TESTCORR')

    def test_set_data(self):
        pass

if __name__ == "__main__":
    unittest.main()
