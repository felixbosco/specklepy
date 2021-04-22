import unittest
import astropy.units as u
from specklepy.mock.target import Target
from specklepy.mock.telescope import Telescope
from specklepy.mock.detector import Detector
from specklepy.scripts.generate import generate_exposure, get_objects


class TestGenerateExposure(unittest.TestCase):

    def setUp(self):
        self.parameterfile = 'specklepy/tests/files/test_synthetic_exposures.par'
        self.star_table = 'specklepy/tests/files/example_star_table_arcsec.dat'
        self.psf_source = 'specklepy/tests/files/psf_short_exposures.fits'
        self.outfile = 'specklepy/tests/files/mock/exposure.fits'

        self.target = Target(band='H', star_table_file=self.star_table, sky_background=14.4)
        self.telescope = Telescope(8.2 * u.m, central_obscuration=0.14, name="VLT Unit Telescope",
                                   psf_source=self.psf_source)
        self.detector = Detector((1024, 1024),
                                 pixel_scale=0.0106 * u.arcsec,
                                 readout_noise=35 * u.electron,
                                 system_gain=17 * u.electron / u.adu,
                                 optics_transmission=0.9,
                                 quantum_efficiency=0.9 * u.electron / u.ph,
                                 saturation_level=7200 * u.electron)
        self.DIT = 0.2 * u.s

    def test_get_objects(self):
        get_objects(self.parameterfile, debug=True)

    def test_generate_exposure(self):
        generate_exposure(self.target, self.telescope, self.detector, self.DIT, number_files=15, n_frames_limit=10,
                          out_file_name=self.outfile, debug=False, dithers=[(10, 0), (-10, 0)], cards={'OBJECT': 'SYNTHETIC', 'OBSTYPE': 'SCIENCE'})


if __name__ == "__main__":
    unittest.main()
