import numpy as np
from os import path
from astropy.io import fits
from astropy.table import Table
# from astropy.nddata import NDData
# from photutils.psf import extract_stars
# from photutils.psf import EPSFBuilder

from holopy.logging import logging
from holopy.io.paramhandler import ParamHandler
from holopy.io.psffile import PSFfile
from holopy.core.aperture import Aperture
from holopy.utils.plot import imshow

class PSFExtraction(object):

    def __init__(self, params):
        if not isinstance(params, ParamHandler):
            raise TypeError("params argument of the PSFExtractor class must be instance of holopy.io.paramhandler.ParamHandler!")
        self.params = params
        self.radius = params.psfRadius

        # Extract stars out of params.refSourceFile
        self.star_table = Table.read(params.refSourceFile, format='ascii')


    @property
    def box_size(self):
        return self.radius * 2 + 1


    def init_ref_apertures(self, filename):
        self.ref_apertures = []
        for star in star_table:
            self.ref_apertures.append(Aperture(star.x, star.y, self.radius, data=filename))

    # def extract(self):
    #
    #
    #     for file in self.params.inFiles:
    #         logging.info("Extracting PSFs from file {}".format(file))
    #         psf_file = PSFfile(file, self.params.tmpDir, frame_shape=(box_size, box_size))
    #
    #         with fits.open(file) as hdulist:
    #             cube = hdulist[0].data
    #
    #             frame_number = hdulist[0].header['NAXIS3']
    #             for frame_index, frame in enumerate(cube):
    #                 print("\rExtracting PSF from frame {}/{}".format(frame_index + 1, frame_number), end='')
    #                 with fits.open(file) as hdulist:
    #
    #
    #                     psf_file.update_frame(frame_index, epsf)
    #             print('\r')
