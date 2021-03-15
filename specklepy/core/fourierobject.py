import numpy as np
from numpy.fft import fft2, ifft2, fftshift
import os
from tqdm import trange

from astropy.io import fits

from specklepy.core.alignment import derive_pad_vectors, pad_array
from specklepy.core.psfmodel import PSFModel
from specklepy.exceptions import SpecklepyValueError
from specklepy.logging import logger
from specklepy.utils.transferfunctions import otf


class FourierObject(object):

    """Reconstruction of the Fourier transformed object.

    This class computes the Fourier transformed object information, as defined in Eq. 1 (Schoedel et al., 2013).
    Therefore it estimates the padding of the PSF and image frames.
    """

    def __init__(self, in_files, psf_files, shifts, mode='same', in_dir=None):
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
        """

        # Assert that there are the same number of inFiles and psfFiles, which should be the case after running the
        # holography function.
        if not len(in_files) == len(psf_files):
            raise ValueError(f"The number of input files ({len(in_files)}) and PSF files ({len(psf_files)}) do not "
                             f"match!")
        self.in_files = in_files
        self.psf_files = psf_files
        self.shifts = shifts

        # Check whether mode is supported
        if mode not in ['same', 'full', 'valid']:
            raise SpecklepyValueError('FourierObject', argname='mode', argvalue=mode,
                                      expected="either 'same', 'full', or 'valid'")
        self.mode = mode
        if in_dir is None:
            self.in_dir = ''
        else:
            self.in_dir = in_dir

        # Extract padding vectors for images and reference image
        logger.info("Initializing padding vectors")
        # files_contain_data_cubes = fits.getdata(in_files[0]).ndim == 3
        self.pad_vectors, self.reference_image_pad_vector = derive_pad_vectors(shifts=shifts,
                                                                               cube_mode=False,
                                                                               return_reference_image_pad_vector=True)
        file_index = 0
        image_pad_vector = self.pad_vectors[file_index]

        # Get example image frame, used as final image size
        image_file = in_files[file_index]
        logger.info(f"\tUsing example image frame from {image_file!r}")
        img = fits.getdata(os.path.join(self.in_dir, image_file))[0]  # Remove time axis padding
        img = pad_array(array=img, pad_vector=image_pad_vector, mode=mode,
                        reference_image_pad_vector=self.reference_image_pad_vector)
        logger.info(f"\tShift: {shifts[file_index]}")
        logger.info(f"\tShape: {img.shape}")

        # Get example PSF frame
        psf_file = psf_files[file_index]
        logger.info(f"\tUsing example PSF frame from {psf_file!r}")
        psf = fits.getdata(psf_file)[0]
        logger.info(f"\tShape: {psf.shape}")

        # Estimate the padding vector for the f_psf frames to have the same xy-extent as f_img
        dx = img.shape[0] - psf.shape[0]
        dy = img.shape[1] - psf.shape[1]
        psf_pad_vector = ((dx // 2, int(np.ceil(dx / 2))), (dy // 2, int(np.ceil(dy / 2))))
        logger.info(f"\tPad_width for PSFs: {psf_pad_vector}")

        # Apply padding to PSF frame
        psf = np.pad(psf, psf_pad_vector, mode='constant', )
        if not img.shape == psf.shape:
            raise ValueError(f"The Fourier transformed images and PSFs have different shape, {img.shape} and "
                             f"{psf.shape}. Something went wrong with the padding!")
        self.psf_pad_vector = psf_pad_vector

        # Initialize the enumerator, denominator and Fourier object attributes
        self.enumerator = np.zeros(img.shape, dtype='complex128')
        self.denominator = np.zeros(img.shape, dtype='complex128')
        self.fourier_image = np.zeros(img.shape, dtype='complex128')

    def coadd_fft(self):
        """Co-add the Fourier transforms of the image and PSF frames.

        Returns:
            fourier_image (np.ndarray, dtype=np.comlex128):
                Fourier-transformed object reconstruction.
        """

        # Padding and Fourier transforming the images
        logger.info("Padding the images and PSFs...")

        for file_index in trange(len(self.in_files), desc="Processing files"):

            # Open PSF and image files
            psf_cube = fits.getdata(self.psf_files[file_index])
            image_cube = fits.getdata(os.path.join(self.in_dir, self.in_files[file_index]))
            n_frames = image_cube.shape[0]

            for frame_index in trange(n_frames, desc="Fourier transforming frames"):

                # Padding and transforming the image
                img = pad_array(array=image_cube[frame_index],
                                pad_vector=self.pad_vectors[file_index],
                                mode=self.mode,
                                reference_image_pad_vector=self.reference_image_pad_vector)
                f_img = fftshift(fft2(img))

                # Padding and Fourier transforming PSF
                psf = psf_cube[frame_index]
                psf = np.pad(psf, self.psf_pad_vector, mode='constant',)
                f_psf = fftshift(fft2(psf))

                # Co-adding for the average
                self.enumerator += np.multiply(f_img, np.conjugate(f_psf))
                self.denominator += np.abs(np.square(f_psf))

        # Compute the object:
        # Note that by this division implicitly does averaging. By this implicit summing up of enumerator and
        # denominator, this computation is cheaper in terms of memory usage
        self.fourier_image = np.divide(self.enumerator, self.denominator)

        from IPython import embed
        embed()

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

        return self.fourier_image

    def ifft(self, total_flux=None):
        """Compute the image by an inverse Fourier transformation of the Fourier-plane image.

        Args:
            total_flux (float, optional):
                Total flux value. The image will be rescaled to obey the total flux, if provided.

        Returns:
            image (np.ndarray):
                Image-plane image.
        """
        logger.info("Inverse Fourier transformation of the object...")
        image = ifft2(self.fourier_image)
        image = np.abs(image)
        if total_flux is not None:
            image_scale = total_flux / np.sum(image)
            image = np.multiply(image, image_scale)
        return image
