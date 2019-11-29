import unittest
import os


class TestApertureAnalysis(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/apertureanalysis.py \
                    -f data/test/example_cube.fits \
                    -i 658 723 \
                    -r 100')

    # def test_first_guess(self):
    #     os.system('python specklepy/scripts/apertureanalysis.py \
    #                 -f data/test/example_cube.fits \
    #                 -i 654 714 \
    #                 -r 100')

    # def test_execute_from_Fourier_file(self):
    #     os.system('python specklepy/scripts/apertureanalysis.py \
    #                 -F data/test/example_cube_Fourier.fits \
    #                 -i 658 723 \
    #                 -r 200')

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
    #         os.system('python specklepy/scripts/apertureanalysis.py \
    #                     -f data/example/esm_200ms_x100_ct.fits \
    #                     -i {} {} \
    #                     -r {} \
    #                     -o data/example/esm_200ms_x100_ct_{}mag_{}pix.dat'.format(*indices[mag], radius, mag, radius))

    # def test_exposure_time(self):
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
    #     mag = '12.0'
    #     exposure_times = [200, 600, 1000, 1400, 1800, 2200]
    #     for expt in exposure_times:
    #         os.system('python specklepy/scripts/apertureanalysis.py \
    #                     -f data/example/glao_{}ms_x100_ct.fits \
    #                     -i {} {} \
    #                     -r {} \
    #                     -o data/example/glao_{}ms_x100_ct_{}mag_{}pix.dat'.format(expt, *indices[mag], radius, expt, mag, radius))

    # def test_aoli_HD(self):
    #     radius = 1.5
    #     fov = 36 # arcsec
    #     detectors_peraxis =  2
    #     pix_per_detector = 1024
    #     pixel_scale = fov / detectors_peraxis / pix_per_detector
    #     radius = int(radius / pixel_scale)
    #     print(radius)
    #
    #     index = (517, 420) # Primary CLOSED
    #     index = (513, 416) # Primary OPEN
    #     # index = (516, 888) # Tertiary
    #     from glob import glob
    #     DATA_PATH = '../../aoli/aoli_fits/'
    #
    #     for path in sorted(glob(DATA_PATH + 'HD207470_OPEN_*.fits')):
    #     # for path in ['HD207470_CLOSED_1_Fourier.fits']:
    #         file = path.split('/')[-1]
    #         os.system('python specklepy/scripts/apertureanalysis.py \
    #                  -f {} \
    #                  -i {} {} \
    #                  -r {} \
    #                  -v True \
    #                  -o data/example/{}.dat'.format(path, *index, radius, file))

    # def test_aoli_HIP(self):
    #     radius = 1.5
    #     fov = 36 # arcsec
    #     detectors_peraxis =  2
    #     pix_per_detector = 1024
    #     pixel_scale = fov / detectors_peraxis / pix_per_detector
    #     radius = int(radius / pixel_scale)
    #     print(radius)
    #
    #     index = (488, 522) # Primary CLOSED
    #     # index = (488, 522) # Primary OPEN
    #     from glob import glob
    #     DATA_PATH = '../../aoli/aoli_fits/'
    #
    #     for file_index, path in enumerate(sorted(glob(DATA_PATH + 'HIP_10644_Flat_mirror_*.fits'))):
    #         visual = '--visual' if file_index==0 else ''
    #         file = path.split('/')[-1]
    #         os.system('python specklepy/scripts/apertureanalysis.py \
    #                  -f {} \
    #                  -i {} {} \
    #                  -r {} \
    #                  {} \
    #                  -o data/example/{}.dat'.format(path, *index, radius, visual, file))


if __name__ == "__main__":
    unittest.main()
