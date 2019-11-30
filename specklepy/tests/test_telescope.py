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

        tel_compute = Telescope(8.2*u.m, psf_source='Gaussian', radius=0.4*u.arcsec, resolution=0.0106*u.arcsec)
        assert np.abs(np.sum(tel_compute.psf) - 1.0) < 1e-6
        if self.visual > 0:
        	imshow(tel_compute.psf)

        tel_airy = Telescope(8.2*u.m, psf_source='AiryDisk', radius=0.1*u.arcsec, resolution=0.0106*u.arcsec)
        if self.visual > 0:
        	imshow(tel_airy.psf)

    def test_call(self):
        pass



if __name__ == "__main__":
    unittest.main()
