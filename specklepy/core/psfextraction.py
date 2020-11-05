import numpy as np
import os
from scipy import ndimage
from tqdm import trange

from astropy.io import fits

from specklepy.core.aperture import Aperture
from specklepy.core.segmentation import Segmentation
from specklepy.logging import logger
from specklepy.io.psffile import PSFFile
from specklepy.io.table import read_table
from specklepy.utils.combine import weighted_mean
from specklepy.plotting.plot import imshow


class ReferenceStars(object):

    """Class that holds a list of reference stars and can extract the PSFs on their positions."""

    def __init__(self, psf_radius, reference_source_file, in_files, save_dir, in_dir=None, field_segmentation=None):
        """

        Args:
            psf_radius (int):
                Radius of the PSF estimate. The PSF frames will have a box size of `2 * psf_radius + 1`, along each
                axes.
            reference_source_file (str):
                Name of the table file that contains the positions of the reference stars.
            in_files (list of str):
                List of input files.
            save_dir (str):
                Directory where the PSF estimates will be stored.
            in_dir (str, optional):
                Path to the input files.
            field_segmentation (list):
                Segmentation of the image into field segments with their own PSF.
        """

        # Store attributes
        self.radius = psf_radius
        self.star_table = read_table(reference_source_file, error=True)
        self.in_files = in_files
        self.save_dir = save_dir
        if in_dir is None:
            self.in_dir = ''
        else:
            self.in_dir = in_dir

        if field_segmentation:
            example_image_shape = fits.getdata(self.in_files[0])[0].shape
            segmentation = Segmentation(*field_segmentation, image_shape=example_image_shape)
            star_positions = []
            for row in self.star_table:
                star_positions.append((row['x'], row['y']))
            if not segmentation.all_covered(positions=star_positions):
                covering = segmentation.all_covered(positions=star_positions, return_all=True)
                raise RuntimeError(f"Image segmentation into ({field_segmentation}) segments requested, but not all "
                                   f"image segments contain reference stars:\n{covering}\nPlease select more reference "
                                   f"or reduce the field segmentation.")
            self.field_segmentation = segmentation

    @property
    def box_size(self):
        return self.radius * 2 + 1

    def init_apertures(self, filename, shift=None):

        if shift is None:
            shift = (0, 0)

        apertures = []
        for star in self.star_table:
            apertures.append(Aperture(star['y'] - shift[0], star['x'] - shift[1], self.radius,
                                      data=os.path.join(self.in_dir, filename), mask='rectangular', crop=True))
        return apertures

    def extract_psfs(self, file_shifts=None, mode='median', align=True, debug=False):
        """Extract the PSF of the list of ReferenceStars frame by frame.

        Long description...

        Args:
            file_shifts (list, optional):
                List of frame shifts for each of the files with respect to the reference file. These will be used to
                adapt the reference star positions. Default is None.
            mode (str, optional):
                Combination mode for PSFs from different apertures.
            align (bool, optional):
                Execute sub-pixel alignments of apertures. Default is True.
            debug (bool, optional):
                Shows the (integrated) apertures if set to True. Default is False.

        Returns:
            psf_files (list):
                List of file names where the PSFs are stored in.
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
        psf_files = []

        # Iterate over input files
        for file_index, file in enumerate(self.in_files):
            # Initialize file by file
            logger.info("Extracting PSFs from file {}".format(file))
            psf_file = PSFFile(file, out_dir=self.save_dir, frame_shape=(self.box_size, self.box_size),
                               in_dir=self.in_dir, header_card_prefix="HIERARCH SPECKLEPY ")
            psf_files.append(psf_file.file_path)

            # Consider alignment of cubes when initializing the apertures, i.e.
            # the position of the aperture in the shifted cube
            if file_shifts is None:
                apertures = self.init_apertures(file)
            else:
                apertures = self.init_apertures(file, shift=file_shifts[file_index])

            # Extract the number of frames in the FITS file from the header
            frame_number = fits.getheader(os.path.join(self.in_dir, file))['NAXIS3']

            # Check apertures visually
            if debug:
                for index, aperture in enumerate(apertures):
                    imshow(aperture.get_integrated(), title="Inspect reference aperture {}".format(index + 1))

            # Extract the PSF by combining the aperture frames in the desired mode
            for frame_index in trange(frame_number, desc="Extracting PSF frame"):
                psfs = np.empty((len(apertures), self.box_size, self.box_size))
                vars = np.ones((len(apertures), self.box_size, self.box_size))
                for aperture_index, aperture in enumerate(apertures):

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

        return psf_files

    def extract_epsfs(self, file_shifts=None, oversampling=4, debug=False):
        """Extract effective PSFs following Anderson & King (2000).

        Args:
            file_shifts (list, optional):
                List of frame shifts for each of the files with respect to the reference file. These will be used to
                adapt the reference star positions. Default is None.
            oversampling (int, optional):
                Factor of oversampling the input pixel grid. Default is 4.
            debug (bool, optional):
                Shows the (integrated) apertures if set to True. Default is False.

        Returns:
            psf_files (list):
                List of file names where the PSFs are stored in.
        """
        # Create a list of psf files and store it to params
        psf_files = []

        # Iterate over input files
        for file_index, file in enumerate(self.in_files):
            # Initialize file by file
            logger.info("Extracting PSFs from file {}".format(file))
            psf_file = PSFFile(file, out_dir=self.save_dir, frame_shape=(self.box_size, self.box_size),
                               in_dir=self.in_dir, header_card_prefix="HIERARCH SPECKLEPY")
            psf_files.append(psf_file.filename)

            # Consider alignment of cubes when initializing the apertures, i.e.
            # the position of the aperture in the shifted cube
            if file_shifts is None:
                apertures = self.init_apertures(file)
            else:
                apertures = self.init_apertures(file, shift=file_shifts[file_index])

            # Extract the number of frames in the FITS file from the header
            frame_number = fits.getheader(os.path.join(self.in_dir, file))['NAXIS3']

            # Extract the PSF by combining the aperture frames in the desired mode
            for frame_index in trange(frame_number, desc="Extracting PSF from frames"):
                
                if debug:
                    if frame_index > 0:
                        break

                # Initialize oversampled grids
                epsf_oversampled = np.zeros((self.box_size * oversampling, self.box_size * oversampling))
                ivar_oversampled = np.zeros((self.box_size * oversampling, self.box_size * oversampling))

                for aperture_index, aperture in enumerate(apertures):
                    xoff = np.floor(aperture.xoffset * oversampling).astype(int) + oversampling // 2
                    yoff = np.floor(aperture.yoffset * oversampling).astype(int) + oversampling // 2

                    # Getting coordinates of aperture and stretching to oversampled image
                    y, x = np.mgrid[0:self.box_size, 0:self.box_size]
                    x *= oversampling
                    y *= oversampling
                    x += xoff
                    y += yoff
                    
                    epsf_oversampled[y, x] += aperture.data[frame_index]
                    ivar_oversampled[y, x] += np.divide(1, aperture.vars)
                
                if debug:
                    imshow(apertures[0].data[frame_index], maximize=True, title=f"Aperture {0}")
                    imshow(epsf_oversampled, maximize=True, title="oversampled ePSF")
                    imshow(ivar_oversampled, maximize=True, title='oversampled IVAR')
                
                # Sample down to the initial grid
                epsf = np.zeros((self.box_size, self.box_size))
                for indizes, value in np.ndenumerate(epsf):
                    weighted_sum = np.multiply(epsf_oversampled, ivar_oversampled)
                    weighted_sum = np.sum(weighted_sum[indizes[0] * oversampling: (indizes[0] + 1) * oversampling,
                                          indizes[1] * oversampling: (indizes[1] + 1) * oversampling])
                    weights_sum = np.sum(ivar_oversampled[indizes[0] * oversampling: (indizes[0] + 1) * oversampling,
                                         indizes[1] * oversampling: (indizes[1] + 1) * oversampling])
                    epsf[indizes] = np.divide(weighted_sum, weights_sum)
                if debug:
                    imshow(epsf, title='ePSF', maximize=True)
        
                psf_file.update_frame(frame_index, epsf)

        return psf_files
