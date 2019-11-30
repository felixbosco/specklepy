import unittest
import astropy.units as u
from specklepy.synthetic.target import Target
from specklepy.utils.plot import imshow


class TestTarget(unittest.TestCase):

    def setUp(self):
        self.star_table = 'data/test/example_star_table_29mas296875.dat'

    def test_init(self):
        Target(band='H', pixel_scale=0.1*u.arcsec, shape=(64, 64), sky_background=13)
        Target(band='H', FoV=(2*u.arcmin, 2*u.arcmin), shape=(64, 64), sky_background=13)
        Target(band='H', pixel_scale=0.1*u.arcsec, FoV=(2*u.arcmin, 2*u.arcmin), star_table=self.star_table)#, sky_background=14)


    def test_call(self):
        target = Target(band='H', pixel_scale=0.1*u.arcsec, FoV=(2*u.arcmin, 2*u.arcmin), star_table=self.star_table, sky_background=13.)
        image = target.get_fieldofview(1024, 0.1*u.arcsec, dither=(1., 0.5))
        imshow(image)



if __name__ == "__main__":
    unittest.main()
