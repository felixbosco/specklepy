import unittest
from holopy.algorithms.matchstarlist import MatchStarList


class TestMatchStarList(unittest.TestCase):

    def setUp(self):
        self.test_files = ['data/test/example_star_list_1.dat',
                            'data/test/example_star_list_2.dat',
                            'data/test/example_star_list_3.dat']

    def test_init(self):
        MatchStarList()

    def test_call(self):
        algorithm = MatchStarList()
        algorithm(self.test_files)


if __name__ == "__main__":
    unittest.main()
