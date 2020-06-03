import unittest
import astropy.units as u
from specklepy.synthetic.target import Target
from specklepy.utils.plot import imshow


class TestTarget(unittest.TestCase):

    def setUp(self):
        self.star_table = 'specklepy/tests/files/example_star_table_29mas296875.dat'

    def test_init(self):
        Target(band='H')
        Target(band='H', sky_background=13)
        Target(band='H', star_table=self.star_table, sky_background=14)


    def test_call(self):
        target = Target(band='H', star_table=self.star_table, sky_background=13.)
        photon_rate_density = target.get_photon_rate_density(field_of_view=30 * u.arcsec, resolution=.5 * u.arcsec, dither=(1., 0.5))
        imshow(photon_rate_density, title='photon_rate_density')



if __name__ == "__main__":
    unittest.main()
