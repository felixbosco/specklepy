import unittest
from holopy.io.parameterset import ParameterSet


class TestParameterSet(unittest.TestCase):

    def setUp(self):
        self.parameter_file = "data/test/test_parfile.ini"
        self.defaults_file = "holopy/config/holography_defaults.cfg"
        self.essential_attributes = ['inDir', 'tmpDir', 'outFile', 'refSourceFile', 'maskRadius', 'noiseThreshold', 'nonsenseKeyWord']
        self.make_dirs = ['inDir', 'tmpDir']


    def test_init(self):
        self.params = ParameterSet(parameter_file=self.parameter_file,
                        defaults_file=self.defaults_file,
                        essential_attributes=self.essential_attributes,
                        make_dirs=self.make_dirs)
        with self.assertRaises(FileNotFoundError):
            ParameterSet(parameter_file=self.parameter_file + "nonsenseExtension", defaults_file=self.defaults_file)


if __name__ == "__main__":
    unittest.main()
