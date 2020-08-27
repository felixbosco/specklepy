from IPython import embed
import numpy as np
from numpy.fft import fft2, ifft2, fftshift
import os
from tqdm import trange

from astropy.io import fits

from specklepy.core.alignment import get_shifts, get_pad_vectors, pad_array
from specklepy.core.aperture import Aperture
from specklepy.core.apodization import apodize
from specklepy.core.psfextraction import ReferenceStars
from specklepy.core.sourceextraction import find_sources
from specklepy.core.ssa import ssa
from specklepy.io.filearchive import FileArchive
from specklepy.io.reconstructionfile import ReconstructionFile
from specklepy.exceptions import SpecklepyValueError
from specklepy.logging import logger
from specklepy.utils.plot import imshow


def holography(params, mode='same', debug=False):
    """Execute the holographic image reconstruction.

    The holographic image reconstruction is an algorithm as outlined, eg. by Schoedel et al (2013, Section 3). This
    function follows that algorithm, see comments in the code. Most of the important functions are imported from other
    modules of specklepy.

    Args:
        params (dict):
            Dictionary that carries all important parameters.
        mode (str, optional):
            Define the size of the output image as 'same' to the reference
            image or expanding to include the 'full' covered field. Default is
            'same'.
        debug (bool, optional):
            Set to True to inspect intermediate results.
            Default is False.

    Returns:
        image (np.ndarray): The image reconstruction.
    """

    logger.info(f"Starting holographic reconstruction...")
    file_archive = FileArchive(file_list=params['PATHS']['inDir'], cards=['DATE'], dtypes=[str])
    in_files = file_archive.files
    in_dir = file_archive.in_dir

    # Input check
    if mode not in ['same', 'full', 'valid']:
        raise SpecklepyValueError('holography()', argname='mode', argvalue=mode,
                                  expected="either 'same', 'full', or 'valid'")

    # Initialize the outfile
    out_file = ReconstructionFile(filename=params['PATHS']['outFile'], files=in_files,
                                  cards={"RECONSTRUCTION": "Holography"}, in_dir=in_dir)

    # (i-ii) Align cubes
    shifts = get_shifts(files=in_files, reference_file=params['PATHS']['alignmentReferenceFile'],
                        lazy_mode=True, return_image_shape=False, in_dir=in_dir, debug=debug)

    # (iii) Compute SSA reconstruction
    image = ssa(in_files, mode=mode, outfile=out_file, in_dir=in_dir, tmp_dir=params['PATHS']['tmpDir'],
                variance_extension_name=params['OPTIONS']['varianceExtensionName'])
    if isinstance(image, tuple):
        # SSA returned a reconstruction image and a variance image
        image, image_var = image
    total_flux = np.sum(image)  # Stored for flux conservation

    # Start iteration from steps (iv) through (xi)
    while True:
        # (iv) Astrometry and photometry, i.e. StarFinder
        find_sources(image=image,
                     fwhm=params['STARFINDER']['starfinderFwhm'],
                     noise_threshold=params['STARFINDER']['noiseThreshold'],
                     background_subtraction=False,
                     writeto=params['PATHS']['allStarsFile'],
                     starfinder='DAO', verbose=False)

        # (v) Select reference stars
        print("\tPlease copy your desired reference stars from the all stars file into the reference star file!")
        input("\tWhen you are done, just hit a key.")

        # (vi) PSF extraction
        ref_stars = ReferenceStars(psf_radius=params['PSFEXTRACTION']['psfRadius'],
                                   reference_source_file=params['PATHS']['refSourceFile'], in_files=in_files,
                                   save_dir=params['PATHS']['tmpDir'], in_dir=in_dir,
                                   field_segmentation=params['PSFEXTRACTION']['fieldSegmentation'])
        if params['PSFEXTRACTION']['mode'].lower() == 'epsf':
            psf_files = ref_stars.extract_epsfs(file_shifts=shifts, debug=debug)
        elif params['PSFEXTRACTION']['mode'].lower() in ['mean', 'median', 'weighted_mean']:
            psf_files = ref_stars.extract_psfs(file_shifts=shifts,
                                                      mode=params['PSFEXTRACTION']['mode'].lower(), debug=debug)
        else:
            raise RuntimeError(f"PSF extraction mode '{params['PSFEXTRACTION']['mode']}' is not understood!")
        logger.info("Saved the extracted PSFs...")

        # (vii) Noise thresholding
        psf_noise_mask = None
        for file in psf_files:
            with fits.open(file, mode='update') as hdu_list:
                n_frames = hdu_list[0].header['NAXIS3']
                if psf_noise_mask is None:
                    psf_noise_mask = get_noise_mask(hdu_list[0].data[0],
                                                    noise_reference_margin=
                                                    params['PSFEXTRACTION']['noiseReferenceMargin'])
                for index in range(n_frames):
                    reference = np.ma.masked_array(hdu_list[0].data[index], mask=psf_noise_mask)
                    background = np.mean(reference)
                    noise = np.std(reference)
                    update = np.maximum(hdu_list[0].data[index] - background -
                                        params['PSFEXTRACTION']['noiseThreshold'] * noise, 0.0)
                    if np.sum(update) == 0.0:
                        raise ValueError("After background subtraction and noise thresholding, no signal is leftover. "
                                         "Please reduce the noiseThreshold!")
                    update = update / np.sum(update)  # Flux sum of order unity
                    hdu_list[0].data[index] = update
                    hdu_list.flush()

        # (viii) Subtraction of secondary sources within the reference apertures
        # TODO: Implement Secondary source subtraction
        pass

        # (ix) Estimate object, following Eq. 1 (Schoedel et al., 2013)
        f_object = FourierObject(in_files, psf_files, shifts=shifts, mode=mode, in_dir=in_dir)
        f_object.coadd_fft()

        # (x) Apodization
        f_object.apodize(params['APODIZATION']['apodizationType'], params['APODIZATION']['apodizationWidth'])

        # (xi) Inverse Fourier transform to retain the reconstructed image
        image = f_object.ifft(total_flux=total_flux)

        # Inspect the latest reconstruction
        if debug:
            imshow(image)

        # Save the latest reconstruction image to outfile
        out_file.data = image

        # Ask the user whether the iteration shall be continued or not
        answer = input("\tDo you want to continue with one more iteration? [yes/no]\n\t")
        if answer.lower() in ['n', 'no']:
            break

    # Repeat astrometry and photometry, i.e. StarFinder on final image
    find_sources(image=image, fwhm=params['STARFINDER']['starfinderFwhm'],
                 noise_threshold=params['STARFINDER']['noiseThreshold'], background_subtraction=False,
                 writeto=params['PATHS']['allStarsFile'], starfinder='DAO', verbose=False)

    # Finally return the image
    return image


def get_noise_mask(frame, noise_reference_margin):
    """Create an annulus-like mask within a given aperture for measuring noise
    and (sky) background.

    Args:
        frame (np.ndarray):
            Image frame within which the mask is derived.
        noise_reference_margin (int):
            Width of the reference annulus in pixels.

    Returns:
        annulus_mask (np.ndarray, dtype=bool):
            Mask array, derived from the frame.
    """
    center = int((frame.shape[0] - 1) / 2)
    radius = center - noise_reference_margin
    tmp = Aperture(center, center, radius, data=frame, crop=False)
    annulus_mask = np.logical_not(tmp.data.mask)
    return annulus_mask


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
        self.pad_vectors, self.reference_image_pad_vector = get_pad_vectors(shifts=shifts,
                                                                            cube_mode=False,
                                                                            return_reference_image_pad_vector=True)
        file_index = 0
        image_pad_vector = self.pad_vectors[file_index]

        # Get example image frame, used as final image size
        image_file = in_files[file_index]
        logger.info(f"\tUsing example image frame from {image_file}")
        img = fits.getdata(os.path.join(self.in_dir, image_file))[0]  # Remove time axis padding
        img = pad_array(array=img, pad_vector=image_pad_vector, mode=mode,
                        reference_image_pad_vector=self.reference_image_pad_vector)
        logger.info(f"\tShift: {shifts[file_index]}")
        logger.info(f"\tShape: {img.shape}")

        # Get example PSF frame
        psf_file = psf_files[file_index]
        logger.info(f"\tUsing example PSF frame from {psf_file}")
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

        return self.fourier_image

    def apodize(self, type, radius):
        """Apodize the Fourier object with a Gaussian kernel.
        
        Args:
            type (str):
                Type of the apodization. Can be either `Gaussian` or `Airy`. See specklepy.core.apodization for details. 
            radius (float):
                Radius of the apodization kernel. This is the standard deviation of a Gaussian kernel or the radius of
                first zero in the case of an Airy function.

        Returns:
            apodized (np.array, dtype=np.complex128):
                Apodized Fourier-plane image.
        """

        logger.info("Apodizing the object...")
        self.fourier_image = apodize(self.fourier_image, type, radius=radius)
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
