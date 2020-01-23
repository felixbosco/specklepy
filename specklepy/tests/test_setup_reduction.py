import unittest
import os


class TestReduce(unittest.TestCase):

    def setUp(self):
        self.path = '"specklepy/tests/files/synthetic/*.fits"'

    def test_call(self):
        os.system('python specklepy/scripts/setup_reduction.py -i specklepy -p {} -o specklepy/tests/files/reduction/files.tab -s FILE'.format(self.path))



if __name__ == "__main__":
    unittest.main()
