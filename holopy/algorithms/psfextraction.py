import numpy as np
from os import path
from astropy.io import fits
from astropy.table import Table
# from astropy.nddata import NDData
# from photutils.psf import extract_stars
# from photutils.psf import EPSFBuilder

from holopy.logging import logging
from holopy.io.paramhandler import ParamHandler
from holopy.io.filehandler import FileHandler
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
        for star in self.star_table:
            # print(star['x'], star['y'], self.radius, filename)
            self.ref_apertures.append(Aperture(star['y'], star['x'], self.radius, data=filename, subset_only=True))


    def extract(self, mode='simple_median'):
        self.params.psfFiles = []
        for file in self.params.inFiles:
            # Initialize file by file
            logging.info("Extracting PSFs from file {}".format(file))
            psf_file = PSFfile(file, outDir=self.params.tmpDir, frame_shape=(self.box_size, self.box_size))
            self.params.psfFiles.append(psf_file.filename)
            self.init_ref_apertures(file)
            frame_number = fits.getheader(file)['NAXIS3']

            # Extract the PSF by combining the aperture frames in the desired mode
            for frame_index in range(frame_number):
                print("\r\tExtracting PSF from frame {}/{}".format(frame_index + 1, frame_number), end='')
                psf = np.empty((len(self.ref_apertures), self.box_size, self.box_size))
                for aperture_index, aperture in enumerate(self.ref_apertures):
                    # Copy aperture into psf
                    psf[aperture_index] = aperture[frame_index]
                    # Normalization of each psf to make median estimate sensible
                    psf[aperture_index] /= np.sum(psf[aperture_index])

                if mode == 'simple_median':
                    psf = np.median(psf, axis=0)
                elif mode == 'aligned_mean':
                    pass
                psf_file.update_frame(frame_index, psf)
            print('\r')
