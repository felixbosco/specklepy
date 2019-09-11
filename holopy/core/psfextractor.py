from os import path
from astropy.table import Table
from astropy.nddata import NDData
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder

from holopy.io.paramhandler import ParamHandler
from holopy.io.psffile import PSFfile

class PSFExtractor(object):

    def __init__(self, params):
        if not isinstance(params, ParamHandler):
            raise TypeError("params argument of the PSFExtractor class must be instance of holopy.io.paramhandler.ParamHandler!")
        self.params = params

        # Extract stars out of params.refSourceFile
        self.star_table = Table.read(params.refSourceFile, format='ascii')


    def extract(self):
        epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, progress_bar=False)

        for file in self.params.inFiles:
            psf_file = PSFfile(file, params.tmpDir, frame_shape=)

        # stars = extract_stars(NDData(data=data), self.star_table, size=25)
        # epsf, fitted_stars = epsf_builder(stars)
