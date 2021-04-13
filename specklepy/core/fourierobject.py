import numpy as np
from numpy.fft import fft2, ifft2, fftshift

from tqdm import trange

from specklepy.core.alignment import FrameAlignment
from specklepy.core.bootstrap import random_draw_vectors
from specklepy.core.psfmodel import PSFModel
from specklepy.exceptions import SpecklepyValueError
from specklepy.io.fits import get_data, get_frame_number
from specklepy.logging import logger
from specklepy.reduction.filter import hot_pixel_mask
from specklepy.utils.array import frame_number
from specklepy.utils.transferfunctions import otf


class FourierObject(object):

    """Reconstruction of the Fourier transformed object.

    This class computes the Fourier transformed object information, as defined in Eq. 1 (Schoedel et al., 2013).
    Therefore it estimates the padding of the PSF and image frames.
    """

    def __init__(self, in_files, psf_files, shifts, mode='same', in_dir=None, custom_mask=None):
        """ Initialize a FourierObject instance.

        Args:
            in_files (list):
                List of paths of the input files.
            psf_files (list):
                List of paths of the PSF files.
            shifts (list):
                List of integer shifts between the files.
            mode (str, optional):
                Define the size of the output image as 'same' to the reference image or expanding to include the 'full'
                covered field. Default is 'same'.
            in_dir (str, optional):
                Path to the input files.
            custom_mask (str or array-like, optional):
                Path to a custom mask or the mask itself, which will be applied to all frames prior to integration.
        """

        # Assert that there are the same number of inFiles and psfFiles, which should be the case after running the
        # holography function.
        if not len(in_files) == len(psf_files):
            raise ValueError(f"The number of input files ({len(in_files)}) and PSF files ({len(psf_files)}) do not "
                             f"match!")
        self.in_files = in_files
        self.psf_files = psf_files
        self.shifts = shifts
        self.in_dir = in_dir

        # Check whether mode is supported
        if mode not in ['same', 'full', 'valid']:
            raise SpecklepyValueError('FourierObject', argname='mode', argvalue=mode,
                                      expected="either 'same', 'full', or 'valid'")
        self.mode = mode

        # Extract padding vectors for images and reference image
        logger.info("Initializing padding vectors")
        self.alignment = FrameAlignment()
        self.alignment.derive_pad_vectors(shifts=shifts)
        file_index = 0

        # Get example image frame, used as final image size
        image_file = in_files[file_index]
        logger.info(f"\tUsing example image frame from {image_file!r}")
        img = get_data(image_file, path=self.in_dir)[0]  # Remove time axis padding
        img = self.alignment.pad_array(array=img, pad_vector_index=file_index, mode=mode)
        logger.info(f"\tShift: {shifts[file_index]}")
        self.shape = img.shape
        logger.info(f"\tShape: {self.shape}")

        # Get example PSF frame
        psf_file = psf_files[file_index]
        logger.info(f"\tUsing example PSF frame from {psf_file!r}")
        psf = get_data(psf_file)[0]
        logger.info(f"\tShape: {psf.shape}")

        # Estimate the padding vector for the f_psf frames to have the same xy-extent as f_img
        dx = self.shape[0] - psf.shape[0]
        dy = self.shape[1] - psf.shape[1]
        psf_pad_vector = ((dx // 2, int(np.ceil(dx / 2))), (dy // 2, int(np.ceil(dy / 2))))
        logger.info(f"\tPad_width for PSFs: {psf_pad_vector}")

        # Apply padding to PSF frame
        psf = np.pad(psf, psf_pad_vector, mode='constant', )
        if not img.shape == psf.shape:
            raise ValueError(f"The Fourier transformed images and PSFs have different shape, {self.shape} and "
                             f"{psf.shape}. Something went wrong with the padding!")
        self.psf_pad_vector = psf_pad_vector

        # Initialize the enumerator, denominator and Fourier object attributes
        self.complex_ratio_image = ComplexRatioImage(self.shape)
        self.fourier_image = None

        # Initialize bootstrap image attribute
        self.bootstrap_images = None
        self.bootstrap_objects = None

        # Set custom mask
        self.custom_mask_file = None
        self._custom_mask = None
        if isinstance(custom_mask, str):
            self.custom_mask_file = custom_mask
        elif isinstance(custom_mask, (bool, np.ndarray)):
            self._custom_mask = custom_mask

    @property
    def number_files(self):
        return len(self.in_files)

    @property
    def custom_mask(self):
        if self._custom_mask is not None:
            return self._custom_mask
        elif self.custom_mask_file is not None:
            return get_data(self.custom_mask_file, dtype=bool)
        else:
            return False

    def coadd_fft(self, fill_value=0, mask_hot_pixels=False, bootstrap=None):
        """Co-add the Fourier transforms of the image and PSF frames.

        Arguments:
            fill_value (float, optional):
                Value to fill masked pixels of the data cube, in case the file contains a mask extension.
            mask_hot_pixels (bool, optional):
                Identify hot pixels in a cube if requested and fill them.
            bootstrap (int, optional):
                Number of bootstrap-resampled images to create in parallel.

        Returns:
            fourier_image (np.ndarray, dtype=np.comlex128):
                Fourier-transformed object reconstruction.
        """

        # Prepare bootstrapping
        if bootstrap is not None:
            self.bootstrap_objects = [ComplexRatioImage(shape=self.shape)] * bootstrap

            # Count total number of frames
            number_frames = 0
            for file in self.in_files:
                number_frames += get_frame_number(file=file, path=self.in_dir)

            # Draw bootstrap coefficients
            bootstrap_draws = random_draw_vectors(number_draws=bootstrap, number_frames=number_frames)
        else:
            bootstrap_draws = None
        frame_counter = 0

        # Padding and Fourier transforming the images
        logger.info("Padding the images and PSFs...")

        # Iterate through files
        for file_index in trange(self.number_files, desc="Processing files"):

            # Open PSF and image files
            psf_cube = get_data(self.psf_files[file_index])
            image_cube = get_data(self.in_files[file_index], path=self.in_dir)
            number_frames = frame_number(image_cube)

            # Extract frame mask
            image_mask = get_data(self.in_files[file_index], path=self.in_dir, extension='MASK', dtype=bool,
                                  ignore_missing_extension=True)
            if image_mask is None:
                image_mask = False

            # Mask hot pixels and combine with image mask
            if mask_hot_pixels:
                # Extract hot pixel mask
                if number_frames == 1:
                    bpm = hot_pixel_mask(image_cube)
                else:
                    bpm = hot_pixel_mask(np.sum(image_cube, axis=0))

                # Combine with image mask
                image_mask = np.logical_or(image_mask, bpm)

            # Combine with the custom mask
            image_mask = np.logical_or(image_mask, self.custom_mask)

            # Iterate through frames
            for frame_index in trange(number_frames, desc="Fourier transforming frames"):

                # Load image and fill masked values with zeros
                frame = image_cube[frame_index]
                frame[image_mask] = fill_value

                # Padding and transforming the image
                img = self.alignment.pad_array(array=frame, pad_vector_index=file_index, mode=self.mode)
                f_img = fftshift(fft2(img))

                # Padding and Fourier transforming PSF
                psf = psf_cube[frame_index]
                psf = np.pad(psf, self.psf_pad_vector, mode='constant',)
                f_psf = fftshift(fft2(psf))

                # Co-adding for the average
                numerator = np.multiply(f_img, np.conjugate(f_psf))
                denominator = np.abs(np.square(f_psf))
                self.complex_ratio_image.numerator += numerator
                self.complex_ratio_image.denominator += denominator

                # Co-add bootstrap images
                if self.bootstrap_objects is not None:
                    for b, bootstrap_object in enumerate(self.bootstrap_objects):
                        bootstrap_coeff = bootstrap_draws[b, frame_counter]
                        if bootstrap_coeff > 0:
                            bootstrap_object.numerator += bootstrap_coeff * numerator
                            bootstrap_object.denominator += bootstrap_coeff * denominator

                # Count frames
                frame_counter += 1

        # Compute the object:
        self.fourier_image = self.complex_ratio_image.evaluate()
        if self.bootstrap_objects is not None:
            self.bootstrap_images = []
            for bootstrap_object in self.bootstrap_objects:
                self.bootstrap_images.append(bootstrap_object.evaluate())

        return self.fourier_image

    def apodize(self, type, radius, crop=False):
        """Apodize the Fourier object with a Gaussian or Airy disk kernel.

        Args:
            type (str):
                Type of the apodization. Can be either `Gaussian` or `Airy`. See specklepy.core.psfmodel for details.
            radius (float):
                Radius of the apodization kernel. This is the standard deviation of a Gaussian kernel or the radius of
                first zero in the case of an Airy function.
            crop (bool, optional):
                Crop corners of the PSF and set them to zero.

        Returns:
            apodized (np.array, dtype=np.complex128):
                Apodized Fourier-plane image.
        """

        # Assert image shape
        if self.fourier_image.shape[0] != self.fourier_image.shape[1]:
            logger.warning("The apodization is applied to a non-quadratic input image. This may cause some "
                           "unpredictable results!")

        logger.info("Apodizing the object...")
        if type is None and radius is None:
            logger.warning(f"Apodization is skipped for either type or radius not being defined!")
            return self.fourier_image

        # Interpret function input and compute apodization PSF
        psf_model = PSFModel(type=type, radius=radius)
        apodization_psf = psf_model(self.fourier_image.shape)

        # Crop corners of the PSF
        if crop:
            threshold = apodization_psf[0, int(psf_model.center[1])]
            apodization_psf -= threshold
            apodization_psf = np.maximum(apodization_psf, 0.0)

        # Normalize to unity
        apodization_psf /= np.sum(apodization_psf)

        # Transform into Fourier space
        apodization_otf = otf(apodization_psf)
        self.fourier_image = np.multiply(self.fourier_image, apodization_otf)

        # Apodize bootstrap images
        if self.bootstrap_images is not None:
            logger.info("Apodizing bootstrap objects...")
            for b, bootstrap_image in enumerate(self.bootstrap_images):
                self.bootstrap_images[b] = np.multiply(bootstrap_image, apodization_otf)

        return self.fourier_image

    def ifft(self, total_flux=None):
        """Compute the image by an inverse Fourier transformation of the Fourier-plane image.

        Args:
            total_flux (float, optional):
                Total flux value. The image will be rescaled to obey the total flux, if provided.

        Returns:
            image (np.ndarray):
                Image-plane image.
            bootstrap_images (list of np.ndarray):
                List of the bootstrap reconstructions of the image. Returns `None`, if no bootstrap objects have been
                computed before.
        """
        logger.info("Inverse Fourier transformation of the object...")
        image = ifft2(self.fourier_image)
        image = np.abs(image)
        if total_flux is not None:
            image_scale = total_flux / np.sum(image)
            image = np.multiply(image, image_scale)

        # Fourier transform of the bootstrap images
        if self.bootstrap_images is not None:
            logger.info("Inverse Fourier transformation of the bootstrap objects...")
            bootstrap_images = []
            for bootstrap_fourier_image in self.bootstrap_images:
                bootstrap_image = ifft2(bootstrap_fourier_image)
                bootstrap_image = np.abs(bootstrap_image)
                if total_flux is not None:
                    bootstrap_image = np.multiply(bootstrap_image, total_flux/np.sum(bootstrap_image))
                bootstrap_images.append(bootstrap_image)
        else:
            bootstrap_images = None

        return image, bootstrap_images


class ComplexRatioImage(object):
    def __init__(self, shape):
        self.shape = shape

        # Initialize numerator and denominator
        self.numerator = np.zeros(shape, dtype='complex128')
        self.denominator = np.zeros(shape, dtype='complex128')

    def evaluate(self):
        """The division implicitly achieves the averaging. By implicitly summing up of enumerator and
        denominator, this computation is cheaper in terms of memory usage."""
        return np.divide(self.numerator, self.denominator)
