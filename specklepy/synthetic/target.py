import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.table import Table

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io import config
from specklepy.io.table import read_table
from specklepy.logging import logger
from specklepy.utils.scaledtuple import ScaledTuple


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
    # Credit https://en.wikipedia.org/wiki/Apparent_magnitude
    # Central wavelength (mum) | Bandpass (mum) |
    photometry_dict = {'U': {'wavelength': 0.36, 'FWHM': 0.15, 'Flux': 1810.0},
                       'B': {'wavelength': 0.44, 'FWHM': 0.22, 'Flux': 4260.0},
                       'V': {'wavelength': 0.55, 'FWHM': 0.16, 'Flux': 3640.0},
                       'R': {'wavelength': 0.64, 'FWHM': 0.23, 'Flux': 3080.0},
                       'I': {'wavelength': 0.79, 'FWHM': 0.19, 'Flux': 2550.0},
                       'J': {'wavelength': 1.26, 'FWHM': 0.16, 'Flux': 1600.0},
                       'H': {'wavelength': 1.6,  'FWHM': 0.23, 'Flux': 1080.0},
                       'K': {'wavelength': 2.22, 'FWHM': 0.23, 'Flux': 670.0},
                       # 'L': {'wavelength': 3.5,  'FWHM': np.nan, 'Flux': np.nan},
                       'g': {'wavelength': 0.52, 'FWHM': 0.14, 'Flux': 3730.0},
                       'r': {'wavelength': 0.67, 'FWHM': 0.14, 'Flux': 4490.0},
                       'i': {'wavelength': 0.79, 'FWHM': 0.16, 'Flux': 4760.0},
                       'z': {'wavelength': 0.91, 'FWHM': 0.13, 'Flux': 4810.0}}

    def __init__(self, band, star_table=None, sky_background=None, photometry_file=None):
        """Instantiate Target class.

        Args:
            band (str):
                Name of the band. Used for extracting the band specific reference flux for magnitude 0.
            star_table (str, optional):
                Name of the file with the data of all stars.
            sky_background (u.Quantity, optional):
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

        if star_table is None or isinstance(star_table, str):
            self.star_table = star_table
            # Read star table already here?
        else:
            raise SpecklepyTypeError('Target', 'star_table', type(star_table), 'str')

        if photometry_file is None or isinstance(photometry_file, str):
            self.photometry_file = photometry_file
        else:
            raise SpecklepyTypeError('Target', 'photometry_file', type(photometry_file), 'str')
        self.band_reference_flux = self.get_reference_flux(self.photometry_file, self.band)

        if isinstance(sky_background, str):
            sky_background = eval(sky_background)
        if sky_background is None:
            self.sky_background_flux = 0.0 / u.Unit('arcsec')**2
        elif isinstance(sky_background, u.Quantity):
            # Interpreting as mag / arcsec**2
            logger.warning("Interpreting sky_background as in units of mag per arcsec**2.")
            self.sky_background_flux = self.magnitude_to_flux(sky_background.value) / u.Unit('arcsec')**2
        elif isinstance(sky_background, (int, float)):
            logger.warning(f"Interpreting scalar type sky_background as {sky_background * u.mag / u.Unit('arcsec')**2}")
            self.sky_background_flux = self.magnitude_to_flux(sky_background) / u.Unit('arcsec')**2
        else:
            raise SpecklepyTypeError('Target', 'sky_background', type(sky_background), 'u.Quantity')

        # Initialize class attributes
        self.shape = None
        self.field_of_view = None
        self.pixel_scale = None
        self.resolution = None
        self.flux_per_pixel = None
        self.stars = None

    @staticmethod
    def from_file(par_file):
        params = config.read(par_file)

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

    def get_reference_flux(self, photometry_file, band, format='ascii'):
        if photometry_file is None:
            fwhm = self.photometry_dict[band]['FWHM']
            flux = self.photometry_dict[band]['Flux'] * u.Unit('Jy')
        else:
            table = Table.read(photometry_file, format=format)
            row_index = np.where(table["Band"] == band)
            fwhm = table['FWHM'][row_index][0]
            flux = table['Flux'][row_index][0] *u.Unit('Jy')
        return (flux / const.h * fwhm * u.Unit('photon')).decompose()

    def magnitude_to_flux(self, magnitude):
        """Convert magnitudes to flux values.

        Args:
            magnitude (int, float, or u.Quantity):
                Magnitude value

        Returns:
            flux (u.Quantity):
                Brightness converted into flux units.
        """
        if isinstance(magnitude, (int, float, np.ndarray)):
            return 10**(magnitude/-2.5) * self.band_reference_flux
        elif isinstance(magnitude, u.Quantity):
            if magnitude.unit != u.Unit('mag'):
                raise SpecklepyValueError('magnitude_to_flux()', 'magnitude unit', magnitude.unit, 'mag')
            else:
                return 10**(magnitude.value/-2.5) * self.band_reference_flux

    def read_star_table(self, file, format='ascii', table_keys=None):
        """Reads a table file and extracts the position and flux of stars.

        Args:
            file (str):
                Name of the file to read in.
            format (str, optional):
                Format of the table file to read. Default is `None`.
            table_keys (dict, optional):
                Keyword dict for 'x' and 'y' position, and 'flux'. Default is `None` and is replaced in the function.

        ToDo:
            * Implemented reading of tables with units.
            * Pass keywords from get_photon_rate_density call to this function.
        """

        # Apply fall bakc values
        if table_keys is None:
            table_keys = {'x': 'x', 'y': 'y', 'flux': 'flux', 'mag': 'mag'}

        # Read table data
        table = read_table(file, format=format)

        # Get data from table columns
        xx = table[table_keys['x']]
        if isinstance(xx, u.Quantity):
            xx = xx.to('arcsec').value

        yy = table[table_keys['y']]
        if isinstance(yy, u.Quantity):
            yy = yy.to('arcsec').value

        if 'mag' in table_keys.keys():
            # Take magnitudes and convert to flux
            magnitudes = table[table_keys['mag']]
            magnitudes = magnitudes.data
            flux = self.magnitude_to_flux(magnitudes)
        else:
            flux = table[table_keys['flux']]

        return Table([xx, yy, flux], names=['x', 'y', 'flux'])

    def get_photon_rate_density(self, field_of_view, resolution, dither=None):
        """Creates an image of the field of view.

        Args:
            field_of_view (u.Quantity or tuple, dtype=u.Quantity):
                Size of the field of view that is covered by the output image.
            resolution (u.Quantity):
                Resolution of the image. Optimally, set it to Telescope.psf_resolution to avoid resampling the image.
            dither (tuple, optional):
                Dither position, relative to the (0, 0) standard phase center.

        Returns:
            photon_rate_density (u.Quantity):
                2D image of the photon rate density towards the standard phase center or dithered position.
        """

        # Input parameters
        if isinstance(field_of_view, u.Quantity):
            self.field_of_view = (field_of_view, field_of_view)
        elif isinstance(field_of_view, tuple):
            self.field_of_view = field_of_view
        elif isinstance(field_of_view, (int, float)):
            logger.warning(f"Interpreting float type FoV as {field_of_view} arcsec")
            field_of_view = field_of_view * u.Unit('arcsec')
            self.field_of_view = (field_of_view, field_of_view)
        else:
            raise SpecklepyTypeError('get_photon_rate_density', 'field_of_view', type(field_of_view), 'tuple')

        # Add 10% FoV to avoid dark margins
        self.field_of_view = (self.field_of_view[0] * 1.1, self.field_of_view[1] * 1.1)

        if isinstance(resolution, (int, float)):
            logger.warning(f"Interpreting float type resolution as {resolution} arcsec")
            resolution = resolution * u.Unit('arcsec')
        elif isinstance(resolution, u.Quantity):
            pass
        else:
            raise SpecklepyTypeError('get_photon_rate_density', 'resolution', type(resolution), 'u.Quantity')
        self.resolution = resolution

        if dither is None:
            phase_center = (0, 0)
        elif isinstance(dither, (tuple, list)):
            if not (isinstance(dither[0], (int, float))):
                raise TypeError("Dithers should be provided as int or float. These are then interpreted as arcsecconds.")
            else:
                phase_center = dither
        else:
            raise SpecklepyTypeError('get_photon_rate_density', 'dither', type(dither), 'tuple')

        # Derive the array shape
        # self.FoV = (self.shape[0] * self.resolution, self.shape[1] * self.resolution)
        shape = (int(self.field_of_view[0] / self.resolution), int(self.field_of_view[1] / self.resolution))
        center = ScaledTuple(shape[0] / 2, shape[1] / 2, scale=self.resolution)
        self.flux_per_pixel = (self.sky_background_flux * self.resolution**2).decompose()

        # Create array with sky background flux
        photon_rate_density = np.ones(shape=shape) * self.flux_per_pixel

        # Add stars from star_table to photon_rate_density
        self.stars = self.read_star_table(self.star_table)
        for row in self.stars:
            # position = (int(center[0] + (row['x'] - phase_center[0]) / self.resolution.to('arcsec').value),
            #             int(center[1] + (row['y'] - phase_center[1]) / self.resolution.to('arcsec').value))
            position = ScaledTuple(row['x'], row['y'], scale=self.resolution.to('arcsec').value, scaled=True,
                                   center=phase_center)
            position.shift(center)
            position = position.index
            flux = row['flux']
            try:
                photon_rate_density.value[position] = np.maximum(photon_rate_density.value[position], flux)
            except IndexError:
                # Star is placed outside the field of view
                pass

        return photon_rate_density
