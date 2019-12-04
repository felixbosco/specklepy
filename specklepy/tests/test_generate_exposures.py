import unittest
import os


class TestGenerateExposure(unittest.TestCase):

    def setUp(self):
        self.parameter_file = 'data/test/test_synthetic_exposures.par'

    def test_call(self):
        os.system('python specklepy/scripts/generate_exposures.py -p {}'.format(self.parameter_file))

    # def test_study_data(self):
    #     path = '../synthetic_observations/par/'
    #     for mode in ['airy']:#, 'sky', 'noao', 'glao']:
    #         for time in ['200', '600', '1200']:
    #             parameter_file = "{}{}_{}ms.par".format(path, mode, time)
    #             os.system('python specklepy/scripts/generate_exposures.py -p {}'.format(parameter_file))



if __name__ == "__main__":
    unittest.main()
