import unittest
import os


class TestScriptHolographicReconstrution(unittest.TestCase):

    def test_execute(self):
        os.system('python holopy/scripts/holographic_reconstruction.py -p {}/data/test/test_parfile.ini'.format(os.getcwd()))


if __name__ == "__main__":
    unittest.main()
