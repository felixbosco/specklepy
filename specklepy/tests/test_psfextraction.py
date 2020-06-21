import unittest
from specklepy.core.psfextraction import ReferenceStars
from specklepy.io.parameterset import ParameterSet


class TestPSFExtraction(unittest.TestCase):

    def setUp(self):
        self.psf_radius = 10
        self.reference_source_file = "specklepy/tests/files/example_ref_sources.dat"
        self.in_files = ["specklepy/tests/files/synthetic/airy_200ms_1.fits",
                         "specklepy/tests/files/synthetic/airy_200ms_2.fits"]
        self.tmp_dir = "specklepy/tests/files/tmp/"
        self.params = {'psf_radius': self.psf_radius, 'reference_source_file': self.reference_source_file,
                       'in_files': self.in_files, 'save_dir': self.tmp_dir}

    def test_init(self):
        ReferenceStars(**self.params)
        ReferenceStars(**self.params, field_segmentation=[2, 1])

    def test_initialize_apertures(self):
        ref_stars = ReferenceStars(**self.params)
        print(ref_stars.star_table)
        ref_stars.init_apertures(self.params['in_files'][0])

    def test_extract_psfs(self):
        ref_stars = ReferenceStars(**self.params)
        ref_stars.extract_psfs()
        ref_stars.extract_psfs(mode='weighted_mean')

    def test_extract_epsfs(self):
        ref_stars = ReferenceStars(**self.params)
        ref_stars.extract_epsfs(debug=False)


if __name__ == "__main__":
    unittest.main()
