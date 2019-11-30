import unittest
import astropy.units as u
from specklepy.synthetic.target import Target
from specklepy.synthetic.telescope import Telescope
from specklepy.synthetic.detector import Detector
from specklepy.synthetic.objectgenerator import get_objects



class TestObjectGenerator(unittest.TestCase):

    def setUp(self):
        self.parameterfile = 'data/test/synthetic_exposures.par'

    def test_call(self):
        get_objects(self.parameterfile, debug=True)



if __name__ == "__main__":
    unittest.main()
