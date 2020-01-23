import unittest
import os


class TestStarfinder(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/starfinder.py \
                    -f specklepy/tests/files/example_cube_ssa.fits \
                    -p specklepy/tests/files/test_reconstruction.par \
                    -o specklepy/tests/files/starfinder_stars.dat')


if __name__ == "__main__":
    unittest.main()
