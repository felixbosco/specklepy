import unittest
import os


class TestReduce(unittest.TestCase):

    def setUp(self):
        self.parameter_file = 'data/test/reduction.par'

    def test_call(self):
        os.system('python specklepy/scripts/reduce.py -p {}'.format(self.parameter_file))



if __name__ == "__main__":
    unittest.main()
