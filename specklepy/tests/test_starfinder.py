import unittest
import os


class TestStarfinder(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/starfinder.py \
                    -f data/test/example_cube_ssa.fits \
                    -p data/test/test_parfile.ini \
                    -o data/test/example_all_sources.dat')


if __name__ == "__main__":
    unittest.main()
