import unittest
import numpy as np
from specklepy.core import holography as holo
from specklepy.io.parameterset import ParameterSet
from specklepy.utils.plot import imshow


class TestHolography(unittest.TestCase):

    def setUp(self):
        self.params = ParameterSet(parameter_file='specklepy/tests/files/test_reconstruction.par',
                        defaults_file='specklepy/config/holography.cfg',
                        essential_attributes={},
                        make_dirs=[])
        self.shifts = [(0, 0)] * 3

    def test_call(self):
        pass

    def test_errors(self):
        with self.assertRaises(ValueError):
            holo.get_fourier_object(self.params, shifts=self.shifts, mode='Nonsense')

if __name__ == "__main__":
    unittest.main()
