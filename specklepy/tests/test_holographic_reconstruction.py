import unittest
import os


class TestHolographicReconstruction(unittest.TestCase):

    # def test_execute(self):
    #     os.system('python specklepy/scripts/holographic_reconstruction.py specklepy/tests/files/test_reconstruction.par')

    def test_execute(self):
        os.system('specklepy holography specklepy/tests/files/test_reconstruction.par')


if __name__ == "__main__":
    unittest.main()
