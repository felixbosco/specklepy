import numpy as np
from scipy import ndimage
from astropy.io import fits
from astropy.table import Table
#from astropy.nddata import NDData
#from astropy.stats import sigma_clipped_stats
#from photutils.psf import extract_stars
#from photutils import EPSFBuilder

from specklepy.logging import logger
from specklepy.io.parameterset import ParameterSet
from specklepy.io.filemanager import FileManager
from specklepy.io.psffile import PSFfile
from specklepy.core.aperture import Aperture
from specklepy.core.combine import weighted_mean
from specklepy.utils.plot import imshow


class ReferenceStars(object):

    """Class that holds a list of reference stars and can extract the PSFs of
    these.

    Long description...
    """

    def __init__(self, params):
        """
        Args:
            params (speckly.io.parameterset.ParameterSet)
        """

        if not isinstance(params, ParameterSet):
            raise TypeError("params argument of the PSFExtractor class must be instance of specklepy.io.parameterset.ParameterSet!")
        self.params = params
        self.radius = params.psfextraction.psfRadius

        # Extract stars out of params.refSourceFile
        self.star_table = Table.read(params.paths.refSourceFile, format='ascii')


    @property
    def box_size(self):
        return self.radius * 2 + 1


    def init_apertures(self, filename, shift=(0, 0)):
        self.apertures = []
        for star in self.star_table:
            self.apertures.append(Aperture(star['y'] - shift[0], star['x'] - shift[1], self.radius, data=filename, mask='rectangular', crop=True, verbose=False))


    def extract_psfs(self, file_shifts=None, mode='median', align=True, debug=False):
        """Extract the PSF of the list of ReferenceStars frame by frame.

        Long description...

        Args:
            mode (str, optional):
            file_shifts (list, optional):
            debug (bool, optional):
                Shows the (integrated) apertures if set to True. Default is
                False.
        """

        # Input parameters
        if mode == 'median':
            func = np.median
        elif mode == 'mean':
            func = np.mean
        elif mode == 'weighted_mean':
            func = weighted_mean
        else:
            raise ValueError('ReferenceStars received unknown mode for extract method ({}).'.format(mode))

        # Create a list of psf files and store it to params
        self.params.psfFiles = []

        # Iterate over params.inFiless
        for file_index, file in enumerate(self.params.inFiles):
            # Initialize file by file
            logger.info("Extracting PSFs from file {}".format(file))
            psf_file = PSFfile(file, outDir=self.params.paths.tmpDir, frame_shape=(self.box_size, self.box_size), header_card_prefix="HIERARCH SPECKLEPY ")
            self.params.psfFiles.append(psf_file.filename)

            # Consider alignment of cubes when initializing the apertures, i.e.
            # the position of the aperture in the shifted cube
            if file_shifts is None:
                file_shift = (0, 0)
            else:
                file_shift = file_shifts[file_index]
            self.init_apertures(file, shift=file_shift)
            frame_number = fits.getheader(file)['NAXIS3']

            # Check apertures visually
            if debug:
                for index, aperture in enumerate(self.apertures):
                    imshow(aperture.get_integrated(), title="Inspect reference aperture {}".format(index + 1))

            # Extract the PSF by combining the aperture frames in the desired mode
            for frame_index in range(frame_number):
                print("\r\tExtracting PSF from frame {}/{}".format(frame_index + 1, frame_number), end='')
                psfs = np.empty((len(self.apertures), self.box_size, self.box_size))
                vars = np.ones((len(self.apertures), self.box_size, self.box_size))
                for aperture_index, aperture in enumerate(self.apertures):

                    flux = aperture[frame_index]
                    var = aperture.vars

                    if align:
                        flux = ndimage.shift(flux, shift=(aperture.xoffset, aperture.yoffset))
                        var = ndimage.shift(var, shift=(aperture.xoffset, aperture.yoffset))

                    # Normalization of each psf to make median estimate sensible
                    psfs[aperture_index] = flux / np.sum(flux)
                    vars[aperture_index] = var / np.sum(flux)

                if mode != 'weighted_mean':
                    psf = func(psfs, axis=0)
                else:
                    psf, var = weighted_mean(psfs, axis=0, vars=vars)

                psf_file.update_frame(frame_index, psf)
            print('\r')


    def extract_epsfs(self, file_shifts=None, oversampling=4, debug=False, **kwargs):
        """Extract effective PSFs following Anderson & King (2000).

        Args:
            file_shifts

        """
        # Create a list of psf files and store it to params
        self.params.psfFiles = []

        # Iterate over params.inFiless
        for file_index, file in enumerate(self.params.inFiles):
            # Initialize file by file
            logger.info("Extracting PSFs from file {}".format(file))
            psf_file = PSFfile(file, outDir=self.params.paths.tmpDir, frame_shape=(self.box_size, self.box_size), header_card_prefix="HIERARCH SPECKLEPY")
            self.params.psfFiles.append(psf_file.filename)

            # Consider alignment of cubes when initializing the apertures, i.e.
            # the position of the aperture in the shifted cube
            if file_shifts is None:
                file_shift = (0, 0)
            else:
                file_shift = file_shifts[file_index]
            self.init_apertures(file, shift=file_shift)

            frame_number = fits.getheader(file)['NAXIS3']

            # Extract the PSF by combining the aperture frames in the desired mode
            for frame_index in range(frame_number):
                print("\r\tExtracting PSF from frame {}/{}".format(frame_index + 1, frame_number), end='')
                
                if debug:
                    if frame_index > 0:
                        break

                # Initialize oversampled grids
                epsf_oversampled = np.zeros((self.box_size * oversampling, self.box_size * oversampling))
                ivar_oversampled = np.zeros((self.box_size * oversampling, self.box_size * oversampling))

                for aperture_index, aperture in enumerate(self.apertures):
                    xoff = np.floor(aperture.xoffset * oversampling).astype(int) + oversampling // 2
                    yoff = np.floor(aperture.yoffset * oversampling).astype(int) + oversampling // 2
                    # print(xoff, yoff)
                    
                    # Getting coordinates of aperture and stretching to oversampled image
                    y, x = np.mgrid[0:self.box_size, 0:self.box_size]
                    x *= oversampling
                    y *= oversampling
                    x += xoff
                    y += yoff
                    
                    epsf_oversampled[y, x] += aperture.data[frame_index]
                    ivar_oversampled[y, x] += np.divide(1, aperture.vars)
                
                if debug:
                    imshow(aperture.data[frame_index], maximize=True, title=f"Aperture {aperture_index}")
                    imshow(epsf_oversampled, maximize=True, title="oversampled ePSF")
                    imshow(ivar_oversampled, maximize=True, title='oversampled IVAR')
                
                # Sample down to the initial grid
                epsf = np.zeros((self.box_size, self.box_size))
                for indizes, value in np.ndenumerate(epsf):
                    weighted_sum = np.multiply(epsf_oversampled, ivar_oversampled)
                    weighted_sum = np.sum(weighted_sum[indizes[0] * oversampling : (indizes[0] + 1) * oversampling , indizes[1] * oversampling : (indizes[1] + 1) * oversampling])
                    weights_sum = np.sum(ivar_oversampled[indizes[0] * oversampling : (indizes[0] + 1) * oversampling , indizes[1] * oversampling : (indizes[1] + 1) * oversampling])
                    epsf[indizes] = np.divide(weighted_sum, weights_sum)
                if debug:
                    imshow(epsf, title='ePSF', maximize=True)
        
                psf_file.update_frame(frame_index, epsf)
            print('\r')


    def extract_epsfs_deprecated(self, file_shifts=None, debug=False):
        """Extract effective PSFs following Anderson & King (2000).

        Args:
            file_shifts

        """
        # Create a list of psf files and store it to params
        self.params.psfFiles = []

        # Iterate over params.inFiless
        for file_index, file in enumerate(self.params.inFiles):

            data = fits.getdata(file)
            if file_shifts is None:
                shift = (0, 0)
            else:
                shift = file_shifts[file_index]
            
            x = self.star_table['x'] - shift[1]
            y = self.star_table['y'] - shift[0]

            mask = ((x > self.radius) & (x < (data.shape[-1] -1 - self.radius)) & (y > self.radius) & (y < (data.shape[-2] -1 - self.radius))) 

            star_table = Table()
            star_table['x'] = x[mask]  
            star_table['y'] = y[mask]  

            print(star_table)

            data = np.sum(data, axis=0)
            mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)  
            data -= median_val 

            nddata = NDData(data=data) 
    
            stars = extract_stars(nddata, star_table, size=self.box_size) 

            if debug:
                for i, star in enumerate(stars):
                    imshow(star.data, title="ePSF star {}".format(i))

            epsf_builder = EPSFBuilder(oversampling=4, maxiters=10, progress_bar=True) 

            epsf, fitted_stars = epsf_builder(stars) 

            if debug:
                imshow(epsf.data, title="ePSF")


