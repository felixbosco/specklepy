import unittest
from holopy.algorithms.aligncubes import AlignCubes


class TestAlignCubes(unittest.TestCase):

    def setUp(self):
        self.test_files = ['data/test/example_shift_image_1.dat',
                            'data/test/example_shift_image_2.dat',
                            'data/test/example_shift_image_3.dat']

    def test_init(self):
        AlignCubes()

    def test_call(self):
        algorithm = AlignCubes()
        algorithm(self.test_files)


if __name__ == "__main__":
    unittest.main()
