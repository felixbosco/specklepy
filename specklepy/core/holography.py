import numpy as np
from numpy.fft import fft2, ifft2, fftshift
from astropy.io import fits

from specklepy.core.alignment import get_shifts, get_pad_vectors, pad_array
from specklepy.core.aperture import Aperture
from specklepy.core.apodization import apodize
from specklepy.core.psfextraction import ReferenceStars
from specklepy.core.sourceextraction import find_sources
from specklepy.core.ssa import ssa
from specklepy.io.parameterset import ParameterSet
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
        params (specklepy.io.parameterset.ParameterSet):
            Class instance that carries all important parameters.
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

    logger.info(f"Starting holographic reconstruction of {len(params.inFiles)} files...")
    if not isinstance(params, ParameterSet):
        logger.warning(f"holography function received params argument of type {type(params)!r} instead of the expected "
                       f"type specklepy.io.parameterset.ParameterSet. This may cause unforeseen errors.")
    if mode not in ['same', 'full', 'valid']:
        raise SpecklepyValueError('holography()', argname='mode', argvalue=mode,
                                  expected="either 'same', 'full', or 'valid'")

    # Initialize the outfile
    params.outFile = ReconstructionFile(filename=params.paths.outFile, files=params.inFiles,
                                        cards={"RECONSTRUCTION": "Holography"})

    # (i-ii) Align cubes
    shifts = get_shifts(files=params.inFiles, reference_file=params.paths.alignmentReferenceFile,
                        lazy_mode=True, return_image_shape=False, debug=debug)

    # (iii) Compute SSA reconstruction
    image = ssa(params.inFiles, mode=mode, outfile=params.outFile, tmp_dir=params.paths.tmpDir,
                variance_extension_name=params.options.varianceExtensionName)
    if isinstance(image, tuple):
        # SSA returned a reconstruction image and a variance image
        image, image_var = image
    total_flux = np.sum(image)  # Stored for flux conservation

    # Start iteration from steps (iv) through (xi)
    while True:
        # (iv) Astrometry and photometry, i.e. StarFinder
        find_sources(image=image,
                     fwhm=params.starfinder.starfinderFwhm,
                     noise_threshold=params.starfinder.noiseThreshold,
                     background_subtraction=False,
                     writeto=params.paths.allStarsFile,
                     starfinder='DAO', verbose=False)

        # (v) Select reference stars
        print("\tPlease copy your desired reference stars from the all stars file into the reference star file!")
        input("\tWhen you are done, just hit a key.")

        # (vi) PSF extraction
        ref_stars = ReferenceStars(psf_radius=params.psfextraction.psfRadius,
                                   reference_source_file=params.paths.refSourceFile, in_files=params.inFiles,
                                   tmp_dir=params.paths.tmpDir)
        if params.psfextraction.mode.lower() == 'epsf':
            params.psf_files = ref_stars.extract_epsfs(file_shifts=shifts, debug=debug)
        elif params.psfextraction.mode.lower() in ['mean', 'median', 'weighted_mean']:
            params.psf_files = ref_stars.extract_psfs(file_shifts=shifts,
                                                      mode=params.psfextraction.mode.lower(), debug=debug)
        logger.info("Saved the extracted PSFs...")

        # (vii) Noise thresholding
        for file in params.psf_files:
            with fits.open(file, mode='update') as hdu_list:
                n_frames = hdu_list[0].header['NAXIS3']
                if not hasattr(params, 'psfNoiseMask'):
                    params.psfNoiseMask = get_noise_mask(hdu_list[0].data[0],
                                                         noise_reference_margin=
                                                         params.psfextraction.noiseReferenceMargin)
                for index in range(n_frames):
                    reference = np.ma.masked_array(hdu_list[0].data[index], mask=params.psfNoiseMask)
                    background = np.mean(reference)
                    noise = np.std(reference)
                    update = np.maximum(hdu_list[0].data[index] - background -
                                        params.psfextraction.noiseThreshold * noise, 0.0)
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
        f_object = FourierObject(params.inFiles, params.psf_files, shifts=shifts, mode=mode)
        f_object.coadd_fft()

        # (x) Apodization
        f_object.apodize(params.apodization.apodizationType, params.apodization.apodizationWidth)

        # (xi) Inverse Fourier transform to retain the reconstructed image
        image = f_object.ifft(total_flux=total_flux)

        # Inspect the latest reconstruction
        if debug:
            imshow(image)

        # Save the latest reconstruction image to outfile
        params.outFile.data = image

        # Ask the user whether the iteration shall be continued or not
        answer = input("\tDo you want to continue with one more iteration? [yes/no]\n\t")
        if answer.lower() in ['n', 'no']:
            break

    # Repeat astrometry and photometry, i.e. StarFinder on final image
    find_sources(image=image, fwhm=params.starfinder.starfinderFwhm, noise_threshold=params.starfinder.noiseThreshold,
                 background_subtraction=False, writeto=params.paths.allStarsFile, starfinder='DAO', verbose=False)

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


def get_fourier_object(in_files, psf_files, shifts, mode='same'):
    """Reconstruction of the Fourier transformed object.

    This function computes the Fourier transformed object information, as defined in Eq. 1 (Schoedel et al., 2013).

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

    Returns:
        f_object (np.ndarray, dtype=complex128):
            Fourier transformed object as a complex128 np.ndarray.
    """

    if mode not in ['same', 'full', 'valid']:
        raise SpecklepyValueError('get_Fourier_object()', argname='mode', argvalue=mode,
                                  expected="either 'same', 'full', or 'valid'")

    files_contain_data_cubes = fits.getdata(in_files[0]).ndim == 3
    pad_vectors, reference_image_pad_vector = get_pad_vectors(shifts=shifts, cube_mode=files_contain_data_cubes,
                                                              return_reference_image_pad_vector=True)

    # Assert that there are the same number of inFiles and psfFiles, which
    # should be the case after running the holography function.
    if not len(in_files) == len(psf_files):
        raise RuntimeError(f"The number of input files ({len(in_files)}) and PSF files ({len(psf_files)}) "
                           f"do not match!")

    # Padding and Fourier transforming the images
    logger.info("Padding the images and PSFs...")
    for file_index, image_file in enumerate(in_files):
        # Initialization
        image_pad_vector = pad_vectors[file_index]
        image_pad_vector.pop(0)
        if file_index == 0:
            # Get final image size
            img = fits.getdata(image_file)[0] # Remove time axis padding
            print(f"\tPadding data from {image_file}")
            img = pad_array(array=img,
                            pad_vector=image_pad_vector,
                            mode=mode,
                            reference_image_pad_vector=reference_image_pad_vector)
            print('\tShift:', shifts[file_index])
            print('\tShape:', img.shape)

            # Get pad vector for PSFs
            psf_file = psf_files[file_index]
            psf = fits.getdata(psf_file)[0]
            # Pad the f_psf cube to have the same xz-extent as f_img
            print(f"\tPadding data from {psf_file}")
            print('\tImage shape:', img.shape)
            print('\tPSF shape:', psf.shape)
            dx = img.shape[0] - psf.shape[0]
            dy = img.shape[1] - psf.shape[1]

            psf_pad_vector = ((dx // 2, int(np.ceil(dx/2))), (dy // 2, int(np.ceil(dy/2))))
            print('\tPad_width:', psf_pad_vector)
            psf = np.pad(psf, psf_pad_vector, mode='constant',)
            try:
                assert img.shape == psf.shape
            except:
                raise ValueError(f"The Fourier transformed images and psfs have different shape, {img.shape} and "
                                 f"{psf.shape}. Something went wrong with the padding!")

            # Initialize the enumerator and denominator
            enumerator = np.zeros(img.shape, dtype='complex128')
            denominator = np.zeros(img.shape, dtype='complex128')

        # Open PSF file
        print(f"\r\tFourier transforming image and PSF file {file_index + 1:4}/{len(in_files):4}", end='')
        psf_cube = fits.getdata(psf_files[file_index])
        for frame_index, frame in enumerate(fits.getdata(image_file)):
            # Padding and transforming the image
            img = pad_array(array=frame,
                            pad_vector=pad_vectors[file_index],
                            mode=mode,
                            reference_image_pad_vector=reference_image_pad_vector)
            f_img = fftshift(fft2(img))

            # Padding and Fourier transforming PSF
            psf = psf_cube[frame_index]
            psf = np.pad(psf, psf_pad_vector, mode='constant',)
            f_psf = fftshift(fft2(psf))

            # Adding for the average
            enumerator += np.multiply(f_img, np.conjugate(f_psf))
            denominator += np.abs(np.square(f_psf))

    print()

    # Compute the object.
    # Note that by this division implicitly does averaging. By this implicit
    # summing up of enumerator and denominator, this computation is cheaper
    # in terms of memory usage
    f_object = np.divide(enumerator, denominator)

    return f_object


class FourierObject(object):

    """Reconstruction of the Fourier transformed object.

    This class computes the Fourier transformed object information, as defined in Eq. 1 (Schoedel et al., 2013).
    Therefore it estimates the padding of the PSF and image frames.
    """

    def __init__(self, in_files, psf_files, shifts, mode='same'):
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
        """

        # Assert that there are the same number of inFiles and psfFiles, which should be the case after running the
        # holography function.
        if not len(in_files) == len(psf_files):
            raise ValueError(f"The number of input files ({len(in_files)}) and PSF files ({len(psf_files)}) do not "\
                               "match!")
        self.in_files = in_files
        self.psf_files = psf_files
        self.shifts = shifts

        # Check whether mode is supported
        if mode not in ['same', 'full', 'valid']:
            raise SpecklepyValueError('FourierObject', argname='mode', argvalue=mode,
                                      expected="either 'same', 'full', or 'valid'")
        self.mode = mode

        # Extract padding vectors for images and reference image
        logger.info("Initializing padding vectors")
        files_contain_data_cubes = fits.getdata(in_files[0]).ndim == 3
        self.pad_vectors, self.reference_image_pad_vector = get_pad_vectors(shifts=shifts,
                                                                            cube_mode=files_contain_data_cubes,
                                                                            return_reference_image_pad_vector=True)
        file_index = 0
        image_pad_vector = self.pad_vectors[file_index]
        image_pad_vector.pop(0)

        # Get example image frame, used as final image size
        image_file = in_files[file_index]
        logger.info(f"\tUsing example image frame from {image_file}")
        img = fits.getdata(image_file)[0]  # Remove time axis padding
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
        # Padding and Fourier transforming the images
        logger.info("Padding the images and PSFs...")
        for file_index, image_file in enumerate(self.in_files):
            # Initialization
            image_pad_vector = self.pad_vectors[file_index]
            image_pad_vector.pop(0)  # Open PSF file

            print(f"\r\tFourier transforming image and PSF file {file_index + 1:4}/{len(self.in_files):4}", end='')
            psf_cube = fits.getdata(self.psf_files[file_index])
            for frame_index, frame in enumerate(fits.getdata(image_file)):
                # Padding and transforming the image
                img = pad_array(array=frame,
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

    def apodize(self, type, width):
        logger.info("Apodizing the object...")
        self.fourier_image = apodize(self.fourier_image, type, radius=width)
        return self.fourier_image

    def ifft(self, total_flux=None):
        logger.info("Inverse Fourier transformation of the object...")
        image = ifft2(self.fourier_image)
        image = np.abs(image)
        if total_flux is not None:
            image_scale = total_flux / np.sum(image)
            image = np.multiply(image, image_scale)
        return image
