import unittest
import os


class TestSSAReconstrution(unittest.TestCase):

    def test_one_cube(self):
        os.system('python specklepy/scripts/ssa_reconstruction.py '
                  '-f specklepy/tests/files/example_cube.fits '
                  '-o specklepy/tests/files/example_cube_ssa.fits')

    def test_multiple_cubes(self):
        os.system('python specklepy/scripts/ssa_reconstruction.py '
                  '-m full '
                  '-f specklepy/tests/files/synthetic/noao_200ms_\*.fits '
                  '-t specklepy/tests/files/tmp/ '
                  '-o specklepy/tests/files/example_cubes_ssa.fits')

    # def test_science_data(self):
    #     os.system('python specklepy/scripts/ssa_reconstruction.py \
    #                 -f "/home/bosco/Documents/sowat/synthetic_observations/reduction/sglao_600ms_?.fits" \
    #                 -t /home/bosco/Documents/sowat/synthetic_observations/tmp/ \
    #                 -o /home/bosco/Documents/sowat/synthetic_observations/ssa_sglao_600ms.fits')

    # def test_airy_data(self):
    #     os.system('python specklepy/scripts/ssa_reconstruction.py \
    #                 -f "/home/bosco/Documents/sowat/synthetic_observations/raw/airy_1200ms_?.fits" \
    #                 -t /home/bosco/Documents/sowat/synthetic_observations/tmp/ \
    #                 -o /home/bosco/Documents/sowat/synthetic_observations/ssa_airy_1200ms.fits')

if __name__ == "__main__":
    unittest.main()
