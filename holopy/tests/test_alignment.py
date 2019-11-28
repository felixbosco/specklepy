import unittest
from holopy.io.filemanager import FileManager
from holopy.core import alignment as align


class TestAlignment(unittest.TestCase):

    def setUp(self):
        self.files = FileManager('/home/bosco/Documents/sowat/simulations/noAO_200ms_x100*_st4_shift.fits')()


    def test_compute_shifts(self):
        align.compute_shifts(self.files)


if __name__ == "__main__":
    unittest.main()
