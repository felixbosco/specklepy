import unittest
from holopy.core.psfextractor import PSFExtractor
from holopy.io.paramhandler import ParamHandler


class TestPSFExtractor(unittest.TestCase):

    def setUp(self):
        self.parameter_file = "data/test/test_parfile.ini"
        self.defaults_file = "holopy/config/holography_defaults.cfg"
        self.essential_attributes = ['inDir', 'tmpDir', 'refSourceFile', 'psfRadius']
        self.make_dirs = ['tmpDir']
        self.params = ParamHandler(parameter_file=self.parameter_file,
                        defaults_file=self.defaults_file,
                        essential_attributes=self.essential_attributes,
                        make_dirs=self.make_dirs)

    def test_init(self):
        PSFExtractor(self.params)

    def test_extract(self):
        psf_extractor = PSFExtractor(self.params)
        psf_extractor.extract()


if __name__ == "__main__":
    unittest.main()
