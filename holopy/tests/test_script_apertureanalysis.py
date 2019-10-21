import unittest
import os


class TestApertureAnalysis(unittest.TestCase):

    def test_execute(self):
        os.system('python holopy/scripts/apertureanalysis.py \
                    -f data/test/example_cube.fits \
                    -i 658 723 \
                    -r 200')

    # def test_first_guess(self):
    #     os.system('python holopy/scripts/apertureanalysis.py \
    #                 -f data/test/example_cube.fits \
    #                 -i 654 714 \
    #                 -r 100')

if __name__ == "__main__":
    unittest.main()
