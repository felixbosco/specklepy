import numpy as np
from scipy import ndimage
from astropy.io import fits
from astropy.table import Table

from specklepy.logging import logging
from specklepy.io.parameterset import ParameterSet
from specklepy.io.filemanager import FileManager
from specklepy.io.psffile import PSFfile
from specklepy.core.aperture import Aperture
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
        self.radius = params.psfRadius

        # Extract stars out of params.refSourceFile
        self.star_table = Table.read(params.refSourceFile, format='ascii')


    @property
    def box_size(self):
        return self.radius * 2 + 1


    def init_apertures(self, filename, shift=(0, 0)):
        self.apertures = []
        for star in self.star_table:
            self.apertures.append(Aperture(star['y'] - shift[0], star['x'] - shift[1], self.radius, data=filename, crop=True, verbose=False))


    def extract_psfs(self, mode='align_median', file_shifts=None, debug=False):
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
        if 'median' in mode:
            combine = np.median
        elif 'mean' in mode:
            combine = np.mean
        else:
            raise ValueError('ReferenceStars received unknown mode for extract method ({}).'.format(mode))

        # Create a list of psf files and store it to params
        self.params.psfFiles = []

        # Iterate over params.inFiless
        for file_index, file in enumerate(self.params.inFiles):
            # Initialize file by file
            logging.info("Extracting PSFs from file {}".format(file))
            psf_file = PSFfile(file, outDir=self.params.tmpDir, frame_shape=(self.box_size, self.box_size), header_prefix="HIERARCH SPECKLEPY ")
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
                psf = np.empty((len(self.apertures), self.box_size, self.box_size))
                for aperture_index, aperture in enumerate(self.apertures):
                    # Copy aperture into psf
                    if 'align' in mode:
                        psf[aperture_index] = ndimage.shift(aperture[frame_index], shift=(aperture.xoffset, aperture.yoffset))
                    elif 'resample' in mode:
                        raise NotImplementedError('Resample mode is not implemented yet!')
                    else:
                        psf[aperture_index] = aperture[frame_index]
                    # Normalization of each psf to make median estimate sensible
                    psf[aperture_index] /= np.sum(psf[aperture_index])

                psf = combine(psf, axis=0)
                psf_file.update_frame(frame_index, psf)
            print('\r')
