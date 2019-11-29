import unittest
import os


class TestScriptSSAReconstrution(unittest.TestCase):

    def test_one_cube(self):
        os.system('python specklepy/scripts/ssa_reconstruction.py \
                    -f data/test/example_cube.fits \
                    -o data/test/example_cube_ssa.fits')

    def test_multiple_cubes(self):
        os.system('python specklepy/scripts/ssa_reconstruction.py \
                    -f "/home/bosco/Documents/sowat/simulations/noAO_200ms_x100*_st4_shift.fits" \
                    -t data/test/tmp/ \
                    -o data/test/example_cubes_ssa.fits')

if __name__ == "__main__":
    unittest.main()
