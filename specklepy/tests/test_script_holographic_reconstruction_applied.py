import unittest
import os


class TestScriptHolographicReconstrution(unittest.TestCase):

    def test_execute(self):
        os.system('python holopy/scripts/holographic_reconstruction.py -p {}/../holography/holography.ini'.format(os.getcwd()))


if __name__ == "__main__":
    unittest.main()
