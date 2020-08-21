import unittest

from specklepy.io.config import read


class TestConfig(unittest.TestCase):

    def setUp(self):
        self.par_file = 'specklepy/tests/files/synthetic/airy_200ms.par'

    def test_read(self):
         config = read(self.par_file)
