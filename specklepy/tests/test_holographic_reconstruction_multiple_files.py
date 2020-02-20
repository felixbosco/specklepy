import unittest
import os


class TestHolographicReconstruction(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/holographic_reconstruction.py specklepy/tests/files/test_holography.par -m full')



if __name__ == "__main__":
    unittest.main()
