from os import path
from astropy.table import Table
from astropy.nddata import NDData
from photutils.psf import extract_stars
from photutils.psf import EPSFBuilder

from holopy.io.paramhandler import ParamHandler

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
            _, filename = path.split(file)
            psf_outfile = filename.replace('.fits', '_psfs.fits')
            psf_outfile = self.params.tmpDir + 'psf/' + psf_outfile
            print(psf_outfile)

        # stars = extract_stars(NDData(data=data), self.star_table, size=25)
        # epsf, fitted_stars = epsf_builder(stars)
