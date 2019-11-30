import unittest
import numpy as np
import astropy.units as u
from specklepy.synthetic.telescope import Telescope
from specklepy.utils.plot import imshow


class TestTelescope(unittest.TestCase):

    def setUp(self):
        self.visual = 1
        self.scao_long_exposure_psf_file = 'data/test/scao_correction_10s_longexposure.fits'
        self.scao_short_exposure_psfs_file = 'data/test/seeing_2ms_shortexposures.fits'

    def test_init(self):
        tel_static = Telescope(8.0*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=self.scao_long_exposure_psf_file )
        assert np.abs(np.sum(tel_static.psf) - 1.0) < 1e-6

        tel_compute = Telescope(8.2*u.m, psf_source='seeing', seeing_fwhm=0.8*u.arcsec, psf_resolution=0.022*u.arcsec, size=128)
        assert np.abs(np.sum(tel_compute.psf) - 1.0) < 1e-6

        if self.visual > 0:
        	imshow(tel_compute.psf)

        # tel_nonstatic = Telescope(8.2*u.m, central_obscuration=0.14, psf_source=self.scao_short_exposure_psfs_file)
        # tel_nonstatic(np.ones((64, 64)), 20*u.mas, integration_time=0.2*u.s)
        # assert np.abs(np.sum(tel_nonstatic.psf) - 1.0) < 1e-6

        tel_airy = Telescope(8.2*u.m, psf_source='airy_model', wavelength=1.63*u.micron)
        # print('Airy model:', tel_airy)

    def test_call(self):
        pass



if __name__ == "__main__":
    unittest.main()
