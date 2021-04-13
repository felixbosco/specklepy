import unittest
import os
from numpy.random import rand

from specklepy.deprecated.psffile import PSFFile


class TestPSFFile(unittest.TestCase):

    def setUp(self):
        self.frame_shape = (25, 25)
        self.path = "specklepy/tests/files/"
        self.file = os.path.join(self.path, 'example_cube.fits')
        self.tmpdir = os.path.join(self.path, 'tmp/')

    def test_init(self):
        PSFFile(in_file=self.file, out_dir=self.tmpdir, frame_shape=self.frame_shape, cards={"RECONSTRUCTION": "Test"}, header_card_prefix="HIERARCH SPECKLEPY ")

    def test_set_data(self):
        psf_file = PSFFile(in_file=self.file, out_dir=self.tmpdir, frame_shape=self.frame_shape, cards={"RECONSTRUCTION": "Test"}, header_card_prefix="HIERARCH SPECKLEPY ")
        psf_file.update_frame(3, rand(*self.frame_shape))


if __name__ == "__main__":
    unittest.main()
