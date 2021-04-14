import numpy as np
import os
from scipy import ndimage
from tqdm import trange
import warnings

from astropy.io import fits

from specklepy.core.aperture import Aperture
from specklepy.core.segmentation import Segmentation
from specklepy.logging import logger
from specklepy.io import FileStream
from specklepy.io.fits import get_frame_number
from specklepy.io.table import read_table
from specklepy.utils.combine import weighted_mean
from specklepy.plotting.utils import imshow


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
        self.in_dir = in_dir if in_dir is not None else ''

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

        # Initialize optional attributes
        self.psf_files = None

    @property
    def box_size(self):
        return self.radius * 2 + 1

    def frame_shape(self, oversampling=1):
        return self.box_size * oversampling, self.box_size * oversampling

    @property
    def number_files(self):
        return len(self.in_files)

    def init_apertures(self, filename, shift=None):

        if shift is None:
            shift = (0, 0)

        apertures = []
        for star in self.star_table:
            apertures.append(Aperture(star['y'] - shift[0], star['x'] - shift[1], self.radius,
                                      file_name=os.path.join(self.in_dir, filename), mask='rectangular', crop=True))
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
            raise ValueError(f"ReferenceStars received unknown mode {mode!r} for extract method!")

        if file_shifts is None:
            file_shifts = [None] * self.number_files

        # Create a list of psf files and store it to params
        self.psf_files = []

        # Iterate over input files
        for file, shift in zip(self.in_files, file_shifts):

            # Extract number of frames
            logger.info(f"Extracting PSFs from file {file!r}")
            frame_number = get_frame_number(file, path=self.in_dir)

            # Initialize PSF file
            psf_file_name = FileStream.default_psf_file_name(file)
            psf_file = FileStream(psf_file_name, path=self.save_dir)
            psf_file.initialize(shape=(frame_number, self.box_size, self.box_size),
                                cards={"SOURCE FILE NAME": os.path.basename(file)})
            self.psf_files.append(psf_file.file_path)

            # Consider alignment of cubes when initializing the apertures, i.e. the position of the aperture in the
            # shifted cube
            apertures = self.init_apertures(file, shift=shift)

            # Check apertures visually
            if debug:
                for index, aperture in enumerate(apertures):
                    imshow(aperture.get_integrated(), title=f"Inspect reference aperture {index + 1}")

            # Extract the PSF by combining the aperture frames in the desired mode
            for frame_index in trange(frame_number, desc="Extracting PSF frame"):
                psfs = np.empty((len(apertures), self.box_size, self.box_size))
                vars = np.ones((len(apertures), self.box_size, self.box_size))
                for aperture_index, aperture in enumerate(apertures):

                    flux = aperture[frame_index]
                    var = aperture.var

                    if align:
                        flux = ndimage.shift(flux, shift=(aperture.y_offset, aperture.x_offset))
                        var = ndimage.shift(var, shift=(aperture.y_offset, aperture.x_offset))

                    # Normalization of each psf to make median estimate sensible
                    psfs[aperture_index] = flux / np.sum(flux)
                    vars[aperture_index] = var / np.sum(flux)

                if mode != 'weighted_mean':
                    psf = func(psfs, axis=0)
                else:
                    psf, var = weighted_mean(psfs, axis=0, vars=vars)

                psf_file.update_frame(frame_index, data=psf)

        return self.psf_files

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

        if file_shifts is None:
            file_shifts = [None] * self.number_files

        # Create a list of psf files and store it to params
        self.psf_files = []

        # Iterate over input files
        for file, shift in zip(self.in_files, file_shifts):

            # Extract frame number
            logger.info(f"Extracting PSFs from file {file!r}")
            frame_number = get_frame_number(file, path=self.in_dir)
            frame_shape = self.frame_shape(oversampling=oversampling)

            # Initialize file by file
            # Initialize PSF file
            psf_file_name = FileStream.default_psf_file_name(file)
            psf_file = FileStream(psf_file_name, path=self.save_dir)
            psf_file.initialize(shape=(frame_number, frame_shape[0], frame_shape[1]),
                                cards={"SOURCE FILE NAME": os.path.basename(file)})
            self.psf_files.append(psf_file.file_path)

            # Consider alignment of cubes when initializing the apertures, i.e. the position of the aperture in the
            # shifted cube
            apertures = self.init_apertures(file, shift=shift)

            # Extract the PSF by combining the aperture frames in the desired mode
            for frame_index in trange(frame_number, desc="Extracting PSF from frames"):
                
                if debug:
                    if frame_index > 0:
                        break

                # Initialize oversampled grids
                epsf_oversampled = np.zeros(frame_shape)
                ivar_oversampled = np.zeros(frame_shape)

                for aperture_index, aperture in enumerate(apertures):
                    xoff = np.floor(aperture.y_offset * oversampling).astype(int) + oversampling // 2
                    yoff = np.floor(aperture.x_offset * oversampling).astype(int) + oversampling // 2

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
        
                psf_file.update_frame(frame_index, data=epsf)

        return self.psf_files

    def create_noise_annulus(self, margin):
        """Create an annulus-like mask within the given aperture for measuring noise and (sky) background.

        Args:
            margin (int):
                Minimum width of the reference annulus in pixels.

        Returns:
            annulus_mask (np.ndarray, dtype=bool):
                Mask array, derived from the frame.
        """

        y, x = np.mgrid[-self.radius: self.radius + 1, -self.radius: self.radius + 1]
        r = np.sqrt(x**2 + y**2)
        return r > self.radius - margin

    def apply_noise_thresholding(self, margin, threshold, extension=0):

        # Create mask for measuring the noise
        noise_reference_region = self.create_noise_annulus(margin=margin)

        # Iterate through PSF files
        for file in self.psf_files:
            psf_file = FileStream(file_name=file)

            # Iterate through frames
            for frame_index in range(get_frame_number(file=file, extension=extension)):

                # Extract statistics from reference annulus
                frame = psf_file.get_frame(frame_index=frame_index, extension=extension)
                reference = frame[noise_reference_region]
                background = np.mean(reference)
                noise = np.std(reference)

                # Subtract background and noise threshold, truncating at zero
                update = np.maximum(frame - background - threshold * noise, 0.0)

                # Normalize
                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    try:
                        update /= np.sum(update)
                    except RuntimeWarning:
                        raise ValueError("After background subtraction and noise thresholding, no signal is leftover. "
                                         "Consider reducing the noiseThreshold!")

                # Update frame in the file
                psf_file.update_frame(frame_index=frame_index, data=update, extension=extension)
