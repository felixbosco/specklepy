import os
import sys
import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.table import Table

from specklepy.logging import logging



class Target(object):

    """ Class carrying information of an astronomical target.

    Long description...
    Note that Target accepts only two out of the three 'shape', 'FoV', and
    'pixel_scale'.

    Attributes:
        shape (tuple, dtype=int):
        FoV (tuple, dtype=Astropy.units.Quantity):
        pixel_scale (astropy.Unit):

    Optional attributes:
        sky_background (int or astropy.units.Quantity, optional):
        config_file (str, optional):
        star_table (str, optional):
        number_stars (int, optional):
    """

    __name__ = 'target'
    typeerror = 'Target received {} argument of {} type, but needs to be {}!'
    # Credit https://en.wikipedia.org/wiki/Photometric_system
    photometry_dict =  {'U': {'wavelength': 0.36, 'FWHM': 0.15, 'Flux': 1810.0},
                        'B': {'wavelength': 0.44, 'FWHM': 0.22, 'Flux': 4260.0},
                        'V': {'wavelength': 0.55, 'FWHM': 0.16, 'Flux': 3640.0},
                        'R': {'wavelength': 0.64, 'FWHM': 0.23, 'Flux': 3080.0},
                        'I': {'wavelength': 0.79, 'FWHM': 0.19, 'Flux': 2550.0},
                        'J': {'wavelength': 1.26, 'FWHM': 0.16, 'Flux': 1600.0},
                        'H': {'wavelength': 1.6,  'FWHM': 0.23, 'Flux': 1080.0},
                        'K': {'wavelength': 2.22, 'FWHM': 0.23, 'Flux': 670.0},
                        'L': {'wavelength': 3.5,  'FWHM': np.nan, 'Flux': np.nan},
                        'g': {'wavelength': 0.52, 'FWHM': 0.14, 'Flux': 3730.0},
                        'r': {'wavelength': 0.67, 'FWHM': 0.14, 'Flux': 4490.0},
                        'i': {'wavelength': 0.79, 'FWHM': 0.16, 'Flux': 4760.0},
                        'z': {'wavelength': 0.91, 'FWHM': 0.13, 'Flux': 4810.0}}



    def __init__(self, band, star_table=None, sky_background=None, photometry_file=None):
        """Instantiate Target class.

        Args:
            band (str): Name of the band. Used for extracting the band specific
                reference flux for magnitude 0.
            star_table (str): Name of the file with the data of all stars.
            sky_background (u.Quantity): Sky background. Int and float inputs
                will be interpreted as mag / arcsec**2.
            photometry_file (str, optional): Name of the file, from which the
                band specific reference flux is extracted.
        """

        # Input parameters
        if isinstance(band, str):
            self.band = band
        else:
            raise TypeError(self.typeerror.format('band', type(band), 'str'))

        if star_table is None:
            self.star_table = None
        elif isinstance(star_table, str):
            self.star_table = star_table
            # Read star table already here?
        else:
            raise TypeError(self.typeerror.format('star_table', type(star_table), 'str'))

        if photometry_file is None or isinstance(photometry_file, str):
            self.photometry_file = photometry_file
        else:
            raise TypeError(self.typeerror.format('photometry_file', type(photometry_file), 'str'))
        self.band_reference_flux = self.get_reference_flux(self.photometry_file, self.band)

        if sky_background is None:
            self.sky_background_flux = 0.0 / u.arcsec**2
        elif isinstance(sky_background, u.Quantity):
    		# Interpreting as mag / arcsec**2
            print("Caution: This function interpretes sky_background as in units of mag per arcsec**2.")
            self.sky_background_flux = self.magnitude_to_flux(sky_background.value) / u.arcsec**2
        elif isinstance(sky_background, int) or isinstance(sky_background, float):
            logging.warning("Interpreting float type sky_background as {}".format(sky_background * u.mag / u.arcsec**2))
            self.sky_background_flux = self.magnitude_to_flux(sky_background) / u.arcsec**2
        else:
            raise TypeError(self.typeerror.format('sky_background', type(sky_background), 'u.Quantity'))



    def __call__(self, *args, **kwargs):
        return self.get_photon_rate_density(*args, **kwargs)



    def __str__(self):
    	tmp = "Target:\n"
    	for key in self.__dict__:
    		if key == 'stars' or key == 'data':
    			continue
    		tmp += "{}: {}\n".format(key, self.__dict__[key])
    	return tmp



    def get_reference_flux(self, photometry_file, band, format='ascii'):
        if photometry_file is None:
            fwhm = self.photometry_dict[band]['FWHM']
            flux = self.photometry_dict[band]['Flux'] * u.Jy
        else:
        	table = Table.read(photometry_file, format=format)
        	row_index = np.where(table["Band"] == band)
        	fwhm = table['FWHM'][row_index][0]
        	flux = table['Flux'][row_index][0] *u.Jy
        return (flux / const.h * fwhm * u.ph).decompose()


    def magnitude_to_flux(self, magnitude):
        """Convert magnitudes to flux values.

        Args:
            magnitude (int, float, or u.Quantity): Magnitude value
        """
        if isinstance(magnitude, int) or isinstance(magnitude, float) or isinstance(magnitude, np.ndarray):
            return 10**(magnitude/-2.5) * self.band_reference_flux
        elif isinstance(magnitude, u.Quantity):
            if magnitude.unit != u.Unit('mag'):
                raise ValueError("The function manitude_to_flux received \
                                magnitude quantity of unit {}, but needs to be \
                                'mag'!".format(magnitude.unit))
            else:
                return 10**(magnitude.value/-2.5) * self.band_reference_flux



    def read_star_table(self, file, format='ascii', keywords=None):
        """Reads a table file and extracts the position and flux of stars.

        Args:
            file (str): Name of the file to read in.
            format (str, optional): Format of the file to read. Passed to
                Table.read(). Default is 'ascii'.
            keywords (dict, optional): Keyword dict for 'x' and 'y' position,
                and 'flux'. Default is None and is replaced in the function.

        To Do:
            * Implemented reading of tables with units.
            * Pass keywords from get_photon_rate_density call to this function.
        """

        # Input parameters
        if keywords is None:
            keywords = {'x': 'x', 'y': 'y', 'flux': 'flux', 'mag': 'mag'}

        table = Table.read(file, format=format)

        # Get data from table columns
        xx = table[keywords['x']]
        yy = table[keywords['y']]
        if 'mag' in keywords.keys():
            # Take magnitudes and convert to flux
            magnitudes = table[keywords['mag']]
            magnitudes = magnitudes.data
            flux = self.magnitude_to_flux(magnitudes)
        else:
            flux = table[keywords['flux']]

        return Table([xx, yy, flux], names=['x', 'y', 'flux'])



    def get_photon_rate_density(self, FoV, resolution, dither=None):
        """Creates an image of the field of view.

        Args:
            FoV (u.Quantity or tuple, dtype=u.Quantity): Size of the field of
                view that is covered by the output image.
            resolution (u.Quantity): Resolution of the image. Optimally, set it
                to Telescope.psf_resolution to avoid resampling the image.
            dither (tuple, optional): Dither position, relative to the (0, 0)
                standard phase center.

        Returns:
            photon_rate_density (u.Quantity): 2D image of the photon rate
                density towards the standard phase center or dithered position.
        """

        # Input parameters
        if isinstance(FoV, u.Quantity):
            self.FoV = (FoV, FoV)
        elif isinstance (FoV, tuple):
            self.FoV = FoV
        elif isinstance(FoV, int) or isinstance(FoV, float):
            logging.warning("Interpreting float type FoV as {}".format(FoV * u.arcsec))
            FoV = FoV * u.arcsec
            self.FoV = (FoV, FoV)
        else:
            return TypeError(self.typeerror.format('shape', type(shape), 'tuple'))
        self.FoV = (self.FoV[0] * 1.1, self.FoV[1] * 1.1) # Add 10% FoV to avoid dark margins

        if isinstance(resolution, int) or isinstance(resolution, float):
            logging.warning("Interpreting float type resolution as {}".format(resolution * u.arcsec))
            resolution = resolution * u.arcsec
        elif isinstance(resolution, u.Quantity):
            pass
        else:
            return TypeError(self.typeerror.format('resolution', type(resolution), 'u.Quantity'))
        self.resolution = resolution

        if dither is None:
            phase_center = (0, 0)
        elif isinstance(dither, tuple):
            phase_center = dither
        else:
            return TypeError(self.typeerror.format('dither', type(dither), 'tuple'))


        # self.FoV = (self.shape[0] * self.resolution, self.shape[1] * self.resolution)
        shape = (int(self.FoV[0] / self.resolution), int(self.FoV[1] / self.resolution))
        center = (shape[0] / 2, shape[1] / 2)
        self.flux_per_pixel = (self.sky_background_flux * self.resolution**2).decompose()

        # Create array with sky background flux
        photon_rate_density = np.ones(shape=shape) * self.flux_per_pixel

        # Add stars from star_table to photon_rate_density
        self.stars = self.read_star_table(self.star_table)
        for row in self.stars:
            position = (int(center[0] + row['x'] - phase_center[0]), int(center[1] + row['y'] - phase_center[1]))
            flux = row['flux']
            try:
                photon_rate_density.value[position] = np.maximum(photon_rate_density.value[position], flux)
            except IndexError:
                # Star is placed outside the field of view
                pass

        return photon_rate_density
