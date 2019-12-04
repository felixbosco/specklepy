import unittest
import os


class TestSSAReconstrution(unittest.TestCase):

    # def test_one_cube(self):
    #     os.system('python specklepy/scripts/ssa_reconstruction.py \
    #                 -f data/test/example_cube.fits \
    #                 -o data/test/example_cube_ssa.fits')
    #
    # def test_multiple_cubes(self):
    #     os.system('python specklepy/scripts/ssa_reconstruction.py \
    #                 -f "/home/bosco/Documents/sowat/simulations/noAO_200ms_x100*_st4_shift.fits" \
    #                 -t data/test/tmp/ \
    #                 -o data/test/example_cubes_ssa.fits')

    def test_science_data(self):
        os.system('python specklepy/scripts/ssa_reconstruction.py \
                    -f "/home/bosco/Documents/sowat/synthetic_observations/reduction/sglao_600ms_?.fits" \
                    -t /home/bosco/Documents/sowat/synthetic_observations/tmp/ \
                    -o /home/bosco/Documents/sowat/synthetic_observations/ssa_sglao_600ms.fits')

    # def test_airy_data(self):
    #     os.system('python specklepy/scripts/ssa_reconstruction.py \
    #                 -f "/home/bosco/Documents/sowat/synthetic_observations/raw/airy_1200ms_?.fits" \
    #                 -t /home/bosco/Documents/sowat/synthetic_observations/tmp/ \
    #                 -o /home/bosco/Documents/sowat/synthetic_observations/ssa_airy_1200ms.fits')

if __name__ == "__main__":
    unittest.main()
