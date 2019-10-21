import unittest
import os


class TestApertureAnalysis(unittest.TestCase):

    # def test_execute(self):
    #     os.system('python holopy/scripts/apertureanalysis.py \
    #                 -f data/test/example_cube.fits \
    #                 -i 658 723 \
    #                 -r 100')

    # def test_first_guess(self):
    #     os.system('python holopy/scripts/apertureanalysis.py \
    #                 -f data/test/example_cube.fits \
    #                 -i 654 714 \
    #                 -r 100')

    def test_execute_from_Fourier_file(self):
        os.system('python holopy/scripts/apertureanalysis.py \
                    -F data/test/example_cube_Fourier.fits \
                    -i 658 723 \
                    -r 200')

    # def test_magnitude(self):
    #     radius = 1.5
    #     pixel_scale = 0.0107421875
    #     radius = int(radius / pixel_scale)
    #
    #     indices = {'12.0': (256, 256),
    #              '12.5': (512, 256),
    #              '13.0': (768, 256),
    #              '13.5': (256, 512),
    #              '14.0': (512, 512),
    #              '14.5': (768, 512),
    #              '15.0': (256, 768),
    #              '15.5': (512, 768),
    #              '16.0': (768, 768)}
    #     for mag in indices:
    #         os.system('python holopy/scripts/apertureanalysis.py \
    #                     -f data/example/esm_200ms_x100_ct.fits \
    #                     -i {} {} \
    #                     -r {} \
    #                     -o data/example/esm_200ms_x100_ct_{}mag_{}pix.dat'.format(*indices[mag], radius, mag, radius))


if __name__ == "__main__":
    unittest.main()
