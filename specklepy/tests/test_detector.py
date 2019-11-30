import unittest
import numpy as np
import astropy.units as u
from specklepy.synthetic.detector import Detector


class TestDetector(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        det = Detector((64, 64), pixel_scale=0.01*u.arcsec)

        det_RON = Detector((64, 64), pixel_scale=0.01*u.arcsec, readout_noise=35*u.electron/u.pix)

        exposure_array = np.ones((64, 64)) * u.ph / u.s
        RO = det_RON(exposure_array, integration_time=200*u.ms, target_FoV=det_RON.FoV)

    def test_call(self):
        pass



if __name__ == "__main__":
    unittest.main()
