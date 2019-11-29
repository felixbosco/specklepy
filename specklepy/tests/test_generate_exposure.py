import unittest
import astropy.units as u
from specklepy.synthetic.target import Target
from specklepy.synthetic.telescope import Telescope
from specklepy.synthetic.detector import Detector
from specklepy.synthetic.generate_exposure import generate_exposure


class TestTarget(unittest.TestCase):

    def setUp(self):
        self.star_table = 'data/test/example_star_table.dat'
        self.psf_source = 'data/test/scao_correction_10s_longexposure.fits'
        # self.psf_source = 'data/test/seeing_2ms_shortexposures.fits'
        self.outdir = 'data/test/exposures/'

        self.target = Target(band='H', FoV=11*u.arcsec*2, shape=(1024, 1024), star_table=self.star_table, sky_background=14.4)
        self.telescope = Telescope(8.2*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=self.psf_source)
        self.detector = Detector((1024, 1024),
        						pixel_size=0.0106*u.arcsec,
        						readout_noise=35*u.electron/u.pix,
        						system_gain=17*u.electron/u.adu,
        						optics_transmission=0.9,
                         		quantum_efficiency=0.9*u.electron/u.ph,
        						saturation_level=7200*u.adu)
        self.DIT = 0.2*u.s

    def test_call(self):
        generate_exposure(self.target, self.telescope, self.detector, self.DIT, number_frames=10, outdir=self.outdir, filename=None, verbose=0, maximum_number_frames_per_file=5, randomkeywordwithoutmeaning=3)



if __name__ == "__main__":
    unittest.main()
