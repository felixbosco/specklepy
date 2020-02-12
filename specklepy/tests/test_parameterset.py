import unittest
from specklepy.io.parameterset import ParameterSet


class TestParameterSet(unittest.TestCase):

    def setUp(self):
        self.parameter_file = "specklepy/tests/files/test_reconstruction.par"
        self.defaults_file = "specklepy/config/holography.cfg"
        self.essential_attributes = {'paths': ['inDir', 'tmpDir', 'outFile', 'refSourceFile'],
                                     'psfextraction': ['psfRadius'],
                                     'starfinder': ['noiseThreshold']}
        self.make_dirs = ['tmpDir']


    def test_init(self):
        self.params = ParameterSet(parameter_file=self.parameter_file,
                        defaults_file=self.defaults_file,
                        essential_attributes=self.essential_attributes,
                        make_dirs=self.make_dirs)
        with self.assertRaises(AttributeError):
            ParameterSet(parameter_file=self.parameter_file, defaults_file=self.defaults_file, essential_attributes={'test': ['nonsenseKeyWord']})
        with self.assertRaises(FileNotFoundError):
            ParameterSet(parameter_file=self.parameter_file + "nonsenseExtension", defaults_file=self.defaults_file)


if __name__ == "__main__":
    unittest.main()
