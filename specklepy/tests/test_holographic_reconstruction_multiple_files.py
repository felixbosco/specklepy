import unittest
import os


class TestHolographicReconstruction(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/holographic_reconstruction.py -p specklepy/tests/files/test_holography.par')



if __name__ == "__main__":
    unittest.main()
