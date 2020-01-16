import unittest
import os


class TestGenerateExposures(unittest.TestCase):

    def setUp(self):
        self.parameter_file = 'data/test/test_synthetic_exposures.par'

    def test_call(self):
        os.system('python specklepy/scripts/generate_exposures.py -p {}'.format(self.parameter_file))

    def test_study_data(self):
        path = 'data/example/synthetic/'
        for mode in ['airy', 'dark', 'sky', 'noao', 'glao']:
            for time in ['200', '600', '1200']:
                parameter_file = "{}{}_{}ms.par".format(path, mode, time)
                os.system('python specklepy/scripts/generate_exposures.py -p {}'.format(parameter_file))
        os.system('python specklepy/scripts/generate_exposures.py -p {}flat.par'.format(path))



if __name__ == "__main__":
    unittest.main()
