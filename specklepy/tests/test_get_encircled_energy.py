import unittest
import os
import glob



class TestGetEncircledEnergy(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/get_encircled_energy.py \
                    -f data/test/example_cube_holo.fits \
                    -i 658 723 \
                    -r 100 \
                    -o data/test/analysis/ \
                    -d True')

    # def test_execute(self):
    #     imagedir = '../synthetic_observations/'
    #     files = glob.glob(imagedir + '*_s*fits')
    #     for file in files:
    #         os.system('python specklepy/scripts/get_encircled_energy.py \
    #                     -f {} \
    #                     -i 660 720 \
    #                     -r 75 \
    #                     -o {}/analysis/ \
    #                     -d True'.format(file, imagedir))




if __name__ == "__main__":
    unittest.main()
