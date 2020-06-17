import unittest
from specklepy.core.reconstruction import Reconstruction


class TestReconstruction(unittest.TestCase):

    def setUp(self):
        self.in_files = []

    def test_init(self):
        Reconstruction(in_files=[], mode='full', outfile='outfile.fits', reference_image=0)
        with self.assertRaises(TypeError):
            Reconstruction(in_files='')
            Reconstruction(in_files=0)
            Reconstruction(in_files=None)
            Reconstruction(in_files=[], mode=None)
            Reconstruction(in_files=[], mode=0)
            Reconstruction(in_files=[], outfile=0)
            Reconstruction(in_files=[], reference_image=0.2)
        with self.assertRaises(ValueError):
            Reconstruction(in_files=[], mode='nonSense')



if __name__ == "__main__":
    unittest.main()
