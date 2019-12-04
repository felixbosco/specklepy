import unittest
import os



class TestApertureAnalysis(unittest.TestCase):

    def test_execute(self):
        os.system('python specklepy/scripts/get_encircled_energy.py \
                    -f data/test/example_cube_holo.fits \
                    -i 658 723 \
                    -r 100 \
                    -o data/test/analysis/ \
                    -d True')



if __name__ == "__main__":
    unittest.main()
