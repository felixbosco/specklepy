import unittest
import os


class TestReduce(unittest.TestCase):

    def setUp(self):
        self.path = '"data/test/synthetic/*.fits"'

    def test_call(self):
        os.system('python specklepy/scripts/setup_reduction.py -i specklepy -p {}'.format(self.path))



if __name__ == "__main__":
    unittest.main()
