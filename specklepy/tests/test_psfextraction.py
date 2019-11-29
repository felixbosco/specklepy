import unittest
from specklepy.algorithms.psfextraction import PSFExtraction
from specklepy.io.parameterset import ParameterSet


class TestPSFExtraction(unittest.TestCase):

    def setUp(self):
        self.parameter_file = "data/test/test_parfile.ini"
        self.defaults_file = "specklepy/config/holography.cfg"
        self.essential_attributes = ['inDir', 'tmpDir', 'refSourceFile', 'psfRadius']
        self.make_dirs = ['tmpDir']
        self.params = ParameterSet(parameter_file=self.parameter_file,
                        defaults_file=self.defaults_file,
                        essential_attributes=self.essential_attributes,
                        make_dirs=self.make_dirs)

    def test_init(self):
        PSFExtraction(self.params)


    def test_initialize_apertures(self):
        algorithm = PSFExtraction(self.params)
        print(algorithm.star_table)
        algorithm.init_ref_apertures(self.params.inFiles[0])


    def test_extract(self):
        algorithm = PSFExtraction(self.params)
        algorithm.extract()
        algorithm.extract(mode='align_median')


if __name__ == "__main__":
    unittest.main()
