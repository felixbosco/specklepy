import unittest
import os
import glob



class TestGetPSFProfile(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/get_psf_profile.py \
                    -f data/test/example_cube_holo.fits \
                    -i 658 723 \
                    -r 100 \
                    -o data/test/analysis/ \
                    -d True')

    # def test_execute_science(self):
    #     imagedir = '../synthetic_observations/'
    #     files = glob.glob(imagedir + '*_s*fits')
    #     for file in files:
    #         os.system('python specklepy/scripts/get_psf_profile.py \
    #                     -f {} \
    #                     -i 660 720 \
    #                     -r 75 \
    #                     -o {}/analysis/ \
    #                     -d True'.format(file, imagedir))

    def test_execute_science_airy(self):
        imagedir = '../synthetic_observations/'
        files = glob.glob(imagedir + 'none_airy_*ms.fits')
        for file in files:
            os.system('python specklepy/scripts/get_psf_profile.py \
                        -f {} \
                        -i 660 720 \
                        -r 75 \
                        -o ../synthetic_observations/analysis/ \
                        -d True'.format(file, imagedir))




if __name__ == "__main__":
    unittest.main()
