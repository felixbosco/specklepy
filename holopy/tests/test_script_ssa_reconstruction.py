import unittest
import os


class TestScriptSSAReconstrution(unittest.TestCase):

    def test_execute(self):
        os.system('python spampy/scripts/extract_spectroastrometric_signal.py \
                    -f data/test/spec2d_cN20170331S0216-pisco_GNIRS_2017Mar31T085412.181.fits \
                    -d Pypeit-Gnirs')


if __name__ == "__main__":
    unittest.main()
