import numpy as np
import os

from astropy.constants import h as PLANCK
from astropy.table import Table
from astropy.units import Unit, Quantity


class PhotometricSystem(object):

    def __init__(self, file_name='default'):

        # Store file name
        if file_name == 'default':
            defaults_file = os.path.join(os.path.dirname(__file__), '../data/photometry.fits')
            self.data_file = os.path.abspath(defaults_file)
        else:
            self.data_file = file_name

        # Initialize attributes
        self._band = None
        self._reference_photon_flux = None

    def reference_photon_flux(self, band=None):
        if (band is None or band == self._band) and self._reference_photon_flux is not None:
            pass
        else:
            self._reference_photon_flux = self.read_reference_photon_flux(band=band)
        return self._reference_photon_flux

    def read_reference_photon_flux(self, band, **kwargs):
        table = Table.read(self.data_file, **kwargs)
        row = table["Band"] == band
        fwhm = table['FWHM'][row][0]
        flux = table['Flux'][row][0]
        return (flux / PLANCK * fwhm * Unit('photon')).decompose()

    def to_photon_flux(self, magnitudes, band=None):
        if isinstance(magnitudes, Quantity):
            if magnitudes.unit == 'mag':
                magnitudes = magnitudes.value

        return np.power(10, np.divide(magnitudes, -2.5)) * self.reference_photon_flux(band=band)
