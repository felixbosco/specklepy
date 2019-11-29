import unittest
from holopy.algorithms.epsfextraction import EPSFExtraction
from holopy.io.parameterset import ParameterSet


class TestEPSFExtraction(unittest.TestCase):

    def setUp(self):
        self.parameter_file = "data/test/test_parfile.ini"
        self.defaults_file = "holopy/config/holography_defaults.cfg"
        self.essential_attributes = ['inDir', 'tmpDir', 'refSourceFile', 'psfRadius']
        self.make_dirs = ['tmpDir']
        self.params = ParameterSet(parameter_file=self.parameter_file,
                        defaults_file=self.defaults_file,
                        essential_attributes=self.essential_attributes,
                        make_dirs=self.make_dirs)

    def test_init(self):
        EPSFExtraction(self.params)

    def test_extract(self):
        epsf_extraction = EPSFExtraction(self.params)
        epsf_extraction.extract()


if __name__ == "__main__":
    unittest.main()
