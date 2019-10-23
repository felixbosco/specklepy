import unittest
import os


class TestScriptSSAReconstrution(unittest.TestCase):

    # def test_execute(self):
        # os.system('python holopy/scripts/ssa_reconstruction.py \
        #             -f data/test/example_cube.fits \
        #             -o data/test/example_cube_ssa.fits')

    def test_aoli(self):
        os.system('python holopy/scripts/ssa_reconstruction.py \
                    -f ../../aoli/aoli_fits/HD207470_CLOSED_1.fits \
                    -o ../../aoli/aoli_fits/HD207470_CLOSED_1_ssa.fits')


if __name__ == "__main__":
    unittest.main()
