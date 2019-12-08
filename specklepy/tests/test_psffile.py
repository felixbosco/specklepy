import unittest
from numpy.random import rand
from specklepy.io.filemanager import FileManager
from specklepy.io.psffile import PSFfile


class TestPSFfile(unittest.TestCase):

    def setUp(self):
        self.frame_shape = (25, 25)

    def test_init(self):
        PSFfile(inFile="data/test/example_cube.fits", outDir="data/test/tmp/", frame_shape=self.frame_shape, cards={"RECONSTRUCTION": "Test"}, header_prefix="HIERARCH SPECKLEPY ")

    def test_set_data(self):
        psf_file = PSFfile(inFile="data/test/example_cube.fits", outDir="data/test/tmp/", frame_shape=self.frame_shape, cards={"RECONSTRUCTION": "Test"}, header_prefix="HIERARCH SPECKLEPY ")
        psf_file.update_frame(3, rand(*self.frame_shape))

if __name__ == "__main__":
    unittest.main()
