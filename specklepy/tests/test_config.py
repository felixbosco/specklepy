import unittest

from specklepy.io import Config


class TestConfig(unittest.TestCase):

    def setUp(self):
        self.par_file = 'specklepy/tests/files/mock/airy_200ms.par'

    def test_read(self):
         config = Config.read(self.par_file)
