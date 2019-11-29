import unittest
import numpy as np
from holopy.core import holography as holo
from holopy.io.parameterset import ParameterSet
from holopy.utils.plot import imshow


class TestHolography(unittest.TestCase):

    def setUp(self):
        self.params = ParameterSet(parameter_file='data/test/test_parfile.ini',
                        defaults_file='holopy/config/holography_defaults.cfg',
                        essential_attributes=[],
                        make_dirs=[])
        self.shifts = [(0, 0)] * 3

    def test_call(self):
        pass

    def test_errors(self):
        with self.assertRaises(ValueError):
            holo.evaluate_object(self.params, shifts=self.shifts, mode='Nonsense')

if __name__ == "__main__":
    unittest.main()
