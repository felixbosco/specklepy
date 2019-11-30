import unittest
import os


class TestGenerateExposure(unittest.TestCase):

    def setUp(self):
        self.parameter_file = 'data/test/synthetic_exposures.par'

    def test_call(self):
        os.system('python specklepy/scripts/generate_exposures.py -p {}'.format(self.parameter_file))



if __name__ == "__main__":
    unittest.main()
