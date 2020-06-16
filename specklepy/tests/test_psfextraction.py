import unittest
from specklepy.core.psfextraction import ReferenceStars
from specklepy.io.parameterset import ParameterSet


class TestPSFExtraction(unittest.TestCase):

    def setUp(self):
        self.parameter_file = "data/test/test_reconstruction.par"
        self.defaults_file = "specklepy/config/holography.cfg"
        self.essential_attributes = ['inDir', 'tmpDir', 'refSourceFile', 'psfRadius']
        self.make_dirs = ['tmpDir']
        self.params = ParameterSet(parameter_file=self.parameter_file, defaults_file=self.defaults_file,
                                   essential_attributes=self.essential_attributes, make_dirs=self.make_dirs)

    def test_init(self):
        ReferenceStars(self.params)

    def test_initialize_apertures(self):
        ref_stars = ReferenceStars(self.params)
        print(ref_stars.star_table)
        ref_stars.init_apertures(self.params.inFiles[0])

    # def test_extract_psfs(self):
    #     ref_stars = ReferenceStars(self.params)
    #     ref_stars.extract_psfs()
    #     # input("Pausing until you want to continue...")
    #     ref_stars.extract_psfs(self.params, mode='weighted_mean')

    def test_extract_epsfs(self):
        ref_stars = ReferenceStars(self.params)
        ref_stars.extract_epsfs(self.params, debug=False)


if __name__ == "__main__":
    unittest.main()
