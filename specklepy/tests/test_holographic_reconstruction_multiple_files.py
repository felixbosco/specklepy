import unittest
import os


class TestHolographicReconstruction(unittest.TestCase):

    def test_execute(self):
        # os.system('python specklepy/scripts/holographic_reconstruction.py -p {}/../holography/holography.ini'.format(os.getcwd()))
        os.system('python specklepy/scripts/holographic_reconstruction.py -p /home/bosco/Documents/sowat/synthetic_observations/holography.par')



if __name__ == "__main__":
    unittest.main()
