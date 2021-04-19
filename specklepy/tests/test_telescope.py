import unittest
import numpy as np
import astropy.units as u
from specklepy.mock.telescope import Telescope
from specklepy.plotting.utils import imshow


class TestTelescope(unittest.TestCase):

    def setUp(self):
        self.visual = 1
        self.scao_long_exposure_psf_file = 'specklepy/tests/files/psf_long_exposure.fits'
        self.scao_short_exposure_psfs_file = 'specklepy/tests/files/psf_short_exposures.fits'
        self.par_file = 'specklepy/tests/files/mock/airy_200ms.par'

    def test_init(self):
        telescope_static_psf = Telescope(8.0*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=self.scao_long_exposure_psf_file )
        assert np.abs(np.sum(telescope_static_psf.psf) - 1.0) < 1e-6

        telescope_gaussian_psf = Telescope(8.2*u.m, psf_source='Gaussian', radius=0.1645*u.arcsec, psf_resolution=0.0106*u.arcsec)
        assert np.abs(np.sum(telescope_gaussian_psf.psf) - 1.0) < 1e-6
        if self.visual > 0:
            imshow(telescope_gaussian_psf.psf, title="Test 'Gaussian' model")

        telescope_airydisk_psf = Telescope(8.2*u.m, psf_source='AiryDisk', radius=0.4777*u.arcsec, psf_resolution=0.0106*u.arcsec)
        if self.visual > 0:
            imshow(telescope_airydisk_psf.psf, title="Test 'AiryDisk' model")

    def test_from_file(self):
        Telescope.from_file(self.par_file)

    # def test_init_static(self):
    #     telescope_static_psf = Telescope(8.0*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=self.scao_long_exposure_psf_file )
    #     assert np.abs(np.sum(telescope_static_psf.psf) - 1.0) < 1e-6
    #
    # def test_init_gaussianmodel(self):
    #     telescope_gaussian_psf = Telescope(8.2*u.m, psf_source='Gaussian', radius=0.1645*u.arcsec, psf_resolution=0.0106*u.arcsec)
    #     assert np.abs(np.sum(telescope_gaussian_psf.psf) - 1.0) < 1e-6
    #     if self.visual > 0:
    #     	imshow(telescope_gaussian_psf.psf, title="Test 'Gaussian' model")
    #
    # def test_init_airymodel(self):
    #     telescope_airydisk_psf = Telescope(8.2*u.m, psf_source='AiryDisk', radius=0.4777*u.arcsec, psf_resolution=0.0106*u.arcsec)
    #     if self.visual > 0:
    #     	imshow(telescope_airydisk_psf.psf, title="Test 'AiryDisk' model")

    def test_call(self):
        pass


if __name__ == "__main__":
    unittest.main()
