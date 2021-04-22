from copy import copy
import numpy as np

from astropy.units import Unit, Quantity
from astropy.table import Column

from specklepy.exceptions import SpecklepyTypeError
from specklepy.io import Config
from specklepy.logging import logger
from specklepy.mock.photometricsystem import PhotometricSystem
from specklepy.mock.startable import StarTable
from specklepy.utils.scaledtuple import Position, ScaledTuple


class Target(object):

    """ Class carrying information of an astronomical target.

    Long description...
    Note that Target accepts only two out of the three 'shape', 'FoV', and
    'pixel_scale'.

    Attributes:
        shape (tuple, dtype=int):
        field_of_view (tuple, dtype=astropy.units.Quantity):
        pixel_scale (astropy.Unit):

    Optional attributes:
        sky_background (int or astropy.units.Quantity, optional):
        config_file (str, optional):
        star_table (str, optional):
        number_stars (int, optional):
    """

    __name__ = 'target'

    def __init__(self, band, star_table_file=None, sky_background=None, photometry_file='default'):
        """Instantiate Target class.

        Args:
            band (str):
                Name of the band. Used for extracting the band specific reference flux for magnitude 0.
            star_table_file (str, optional):
                Name of the file with the data of all stars.
            sky_background (Quantity, optional):
                Sky background. Int and float inputs will be interpreted as mag / arcsec**2.
            photometry_file (str, optional):
                Name of the file, from which the band specific reference flux is extracted.
        """

        # Input parameters
        if isinstance(band, str):
            try:
                self.band = eval(band)
            except NameError:
                self.band = band
        else:
            raise SpecklepyTypeError('Target', 'band', type(band), 'str')

        self.photometric_system = PhotometricSystem(photometry_file)

        if isinstance(sky_background, str):
            sky_background = eval(sky_background)
        if sky_background is None:
            self.sky_background_flux = 0.0 / Unit('arcsec')**2
        elif isinstance(sky_background, Quantity):
            # Interpreting as mag / arcsec**2
            logger.warning("Interpreting sky_background as in units of mag per arcsec**2.")
            self.sky_background_flux = self.photometric_system.to_photon_flux(sky_background.value, band=self.band) / \
                                       Unit('arcsec')**2
        elif isinstance(sky_background, (int, float)):
            logger.warning(f"Interpreting scalar type sky_background as {sky_background * Unit('mag / arcsec**2')}")
            self.sky_background_flux = self.photometric_system.to_photon_flux(sky_background, band=self.band) / \
                                       Unit('arcsec**2')
        else:
            raise SpecklepyTypeError('Target', 'sky_background', type(sky_background), 'Quantity')

        # Initialize class attributes
        self.shape = None
        self.field_of_view = None
        self.pixel_scale = None
        self.resolution = None
        self.flux_per_pixel = None

        self.star_table = None
        if star_table_file is not None:
            self.star_table = StarTable.from_file(star_table_file)
        # self.star_table = self.read_star_table(self.star_table_file)

    @staticmethod
    def from_file(par_file):
        params = Config.read(par_file).params

        # Try known first-level keys
        for key in ['TARGET', 'Target', 'target']:
            if key in params.keys():
                return Target(**params[key])

        # Try full parameter set
        try:
            return Target(**params)
        except TypeError:
            raise RuntimeError(f"Could not identify parameters for initializing a Target instance in parameter "
                               f"file {par_file}")

    def __call__(self, *args, **kwargs):
        return self.get_photon_rate_density(*args, **kwargs)

    def __str__(self):
        tmp = "Target:"
        for key in self.__dict__:
            if key == 'stars' or key == 'data':
                continue
            tmp += f"\n{key}: {self.__dict__[key]}"
        return tmp

    def get_photon_rate_density(self, field_of_view, resolution, dither=None):
        """Creates an image of the field of view.

        Args:
            field_of_view (Quantity or tuple, dtype=Quantity):
                Size of the field of view that is covered by the output image.
            resolution (Quantity):
                Resolution of the image. Optimally, set it to Telescope.psf_resolution to avoid resampling the image.
            dither (tuple, optional):
                Dither position, relative to the (0, 0) standard phase center.

        Returns:
            photon_rate_density (Quantity):
                2D image of the photon rate density towards the standard phase center or dithered position.
        """

        # Input parameters
        if isinstance(resolution, (int, float)):
            logger.warning(f"Interpreting float type resolution as {resolution} arcsec")
            resolution = resolution * Unit('arcsec')
        elif isinstance(resolution, Quantity):
            pass
        else:
            raise SpecklepyTypeError('get_photon_rate_density', 'resolution', type(resolution), 'Quantity')
        self.resolution = resolution

        if isinstance(field_of_view, (int, float)):
            logger.warning(f"Interpreting float type FoV as {field_of_view} arcsec")
            field_of_view = field_of_view * Unit('arcsec')
        elif isinstance(field_of_view, (tuple, list, Quantity)):
            pass
        else:
            raise SpecklepyTypeError('get_photon_rate_density', 'field_of_view', type(field_of_view), 'tuple')

        # Add 10% FoV to avoid dark margins
        self.field_of_view = ScaledTuple(field_of_view, scale=resolution, scaled=True)
        self.field_of_view *= 1.1

        if dither is None:
            phase_center = (0, 0)
        elif isinstance(dither, (tuple, list)):
            if not (isinstance(dither[0], (int, float))):
                raise TypeError("Dithers should be provided as int or float. These are then interpreted as arcsec.")
            else:
                phase_center = dither
        else:
            raise SpecklepyTypeError('get_photon_rate_density', 'dither', type(dither), 'tuple')

        # Define image center for centering star positions around the image center
        center = copy(self.field_of_view) / 2

        # Create array with sky background flux
        self.flux_per_pixel = (self.sky_background_flux * self.resolution**2).decompose()
        photon_rate_density = np.ones(shape=self.field_of_view.index) * self.flux_per_pixel

        # Add stars from star_table to photon_rate_density
        if self.star_table.table is not None:

            distance = 10 * Unit('kpc')
            star_table = self.star_table.observe_at_distance(distance=distance, filter_band=self.band)
            flux_values = self.photometric_system.to_photon_flux(star_table[self.band], band=self.band)
            star_table.add_column(Column(name='flux', data=flux_values))

            for row in star_table:
                position = Position(row['x'], row['y'], scale=self.resolution.to('arcsec').value, scaled=False)
                position.offset(center)
                position.offset(phase_center, scaled=True)
                flux = row['flux']
                try:
                    # photon_rate_density.value[position.index] = np.maximum(photon_rate_density.value[position.index],
                    # flux)
                    photon_rate_density.value[position.index] = photon_rate_density.value[position.index] + flux
                except IndexError:
                    # Star is placed outside the field of view
                    pass

        return photon_rate_density
