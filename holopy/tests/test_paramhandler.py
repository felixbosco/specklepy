import unittest
from holopy.io.paramhandler import ParamHandler


class TestParamHandler(unittest.TestCase):

    def setUp(self):
        self.parameter_file = "data/test/test_parfile.ini"
        self.defaults_file = "holopy/config/holography_defaults.cfg"

    def test_init(self):
        ParamHandler(parameter_file=self.parameter_file, defaults_file=self.defaults_file)

if __name__ == "__main__":
    unittest.main()
