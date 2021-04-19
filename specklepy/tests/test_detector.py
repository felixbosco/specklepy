import unittest
import numpy as np
import astropy.units as u
from specklepy.mock.detector import Detector


class TestDetector(unittest.TestCase):

    def setUp(self):
        self.photon_rate_array = np.ones((128, 128)) * u.ph / u.s
        self.par_file = 'specklepy/tests/files/mock/airy_200ms.par'

    def test_init(self):
        det = Detector((64, 64), pixel_scale=0.01*u.arcsec)
        det_RON = Detector((64, 64), pixel_scale=0.01*u.arcsec, readout_noise=35*u.electron)

    def test_call(self):
        det = Detector((64, 64), pixel_scale=0.01*u.arcsec, readout_noise=35*u.electron)
        det(self.photon_rate_array, integration_time=200*u.ms, photon_rate_resolution=det.pixel_scale)
        det(self.photon_rate_array, integration_time=200*u.ms, photon_rate_resolution=det.pixel_scale / 2)

    def test_from_file(self):
        Detector.from_file(self.par_file)


if __name__ == "__main__":
    unittest.main()
