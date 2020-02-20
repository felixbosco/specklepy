import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
import os
import glob
from datetime import datetime
from astropy.io import fits
from IPython import embed

from specklepy.core.alignment import get_shifts, get_pad_vectors, pad_array
from specklepy.core.aperture import Aperture
from specklepy.core.apodization import apodize
from specklepy.core.psfextraction import ReferenceStars
from specklepy.core.sourceextraction import find_sources
from specklepy.core.ssa import ssa
from specklepy.io.outfile import Outfile
from specklepy.io.parameterset import ParameterSet
from specklepy.io.reconstructionfile import ReconstructionFile
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.logging import logger
from specklepy.utils.plot import imshow
from specklepy.utils.transferfunctions import otf



def holography(params, mode='same', debug=False):
    """Execute the holographic image reconstruction.

    The holographic image reconstruction is an algorithm as outlined, eg. by
    Schoedel et al (2013, Section 3). This function follows that algorithm, see
    comments in the code. Most of the important functions are imported from
    other modules of specklepy.

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

    logger.info("Starting holographic reconstruction of {} files...".format(len(params.inFiles)))
    if not isinstance(params, ParameterSet):
        logger.warning(f"holography function received params argument of type {type(params)!r} instead of the expected type " \
                        "specklepy.io.parameterset.ParameterSet. This may cause unforeseen errors.")
    if mode not in ['same', 'full', 'valid']:
        raise SpecklepyValueError('holography()', argname='mode', argvalue=mode, expected="either 'same', 'full', or 'valid'")


    # Initialize the outfile
    params.outFile = ReconstructionFile(filename=params.paths.outFile, files=params.inFiles, cards={"RECONSTRUCTION": "Holography"})

    # (i-ii) Align cubes
    shifts = get_shifts(files=params.inFiles,
                            reference_file=params.paths.alignmentReferenceFile,
                            lazy_mode=True,
                            return_image_shape=False,
                            debug=debug)

    # (iii) Compute SSA reconstruction
    image = ssa(params.inFiles, mode=mode, outfile=params.outFile, tmp_dir=params.paths.tmpDir, variance_extension_name=params.uncertainties.varianceExtensionName)
    if isinstance(image, tuple):
        # SSA returned a reconstruction image and a variance image
        image, image_var = image
    total_flux = np.sum(image) # Stored for flux conservation

    # Start iteration from steps (iv) through (xi)
    while True:
        # (iv) Astrometry and photometry, i.e. StarFinder
        find_sources(image=image, fwhm=params.starfinder.starfinderFwhm, noise_threshold=params.starfinder.noiseThreshold,
            background_subtraction=False, writeto=params.paths.allStarsFile, starfinder='DAO', verbose=False)

        # (v) Select reference stars
        print("\tPlease copy your desired reference stars from the all stars file into the reference star file!")
        input("\tWhen you are done, just hit a key.")

        # (vi) PSF extraction
        refStars = ReferenceStars(params)
        if params.psfextraction.mode.lower() == 'epsf':
            refStars.extract_epsfs(file_shifts=shifts, debug=debug)
        elif params.psfextraction.mode.lower() in ['mean', 'median', 'weighted_mean']:
            refStars.extract_psfs(file_shifts=shifts, mode=params.psfextraction.mode.lower(), debug=debug)
        logger.info("Saved the extracted PSFs...")

        # (vii) Noise thresholding
        for file in params.psfFiles:
            with fits.open(file, mode='update') as hdulist:
                numberFrames = hdulist[0].header['NAXIS3']
                if not hasattr(params, 'psfNoiseMask'):
                    params.psfNoiseMask = get_noise_mask(hdulist[0].data[0], noise_reference_margin=params.psfextraction.noiseReferenceMargin)
                for index in range(numberFrames):
                    reference = np.ma.masked_array(hdulist[0].data[index], mask=params.psfNoiseMask)
                    background = np.mean(reference)
                    noise = np.std(reference)
                    update = np.maximum(hdulist[0].data[index] - background - params.psfextraction.noiseThreshold * noise, 0.0)
                    if np.sum(update) == 0.0:
                        raise ValueError("After background subtraction and noise thresholding, no signal is leftover. Please reduce the noiseThreshold!")
                    update = update / np.sum(update) # Flux sum of order unity
                    hdulist[0].data[index] = update
                    hdulist.flush()

        # (viii) Subtraction of secondary sources within the reference apertures
        # THIS IS NOT IMPLEMENTED YET!
        pass

        # (ix) Estimate object, following Eq. 1 (Schoedel et al., 2013)
        Fobject = get_Fourier_object(params, shifts=shifts, mode=mode)

        # (x) Apodization
        logger.info("Apodizing the object...")
        Fobject = apodize(Fobject, params.apodization.apodizationType, radius=params.apodization.apodizationWidth)

        # (xi) Inverse Fourier transform to retain the reconstructed image
        image = ifft2(Fobject)
        image = np.abs(image)
        image_scale = total_flux / np.sum(image)
        image = np.multiply(image, image_scale) # assert flux conservation


        # Inspect the latest reconstrction
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



def get_Fourier_object(params, shifts, mode='same'):
    """Reconstruction of the Fourier transformed object.

    This function computes the Fourier transformed object information, as
    defined in Eq. 1 (Schoedel et al., 2013).

    Args:
        params (specklepy.io.parameterset.ParameterSet):
            Class instance that carries all important parameters.
        shifts (list):
            List of integer shifts between the files.
        mode (str, optional):
            Define the size of the output image as 'same' to the reference
            image or expanding to include the 'full' covered field. Default is
            'same'.

    Returns:
        Fobject (np.ndarray, dtype=complex128): Fourier transformed object as a
            complex128 np.ndarray.
    """

    if not isinstance(params, ParameterSet):
        logger.warn("specklepy.core.holography.get_Fourier_object received params argument \
                        of type <{}> instead of the expected type \
                        specklepy.io.parameterset.ParameterSet. This may cause \
                        unforeseen errors.".format(type(params)))
    if mode not in ['same', 'full', 'valid']:
        raise SpecklepyValueError('get_Fourier_object()', argname='mode', argvalue=mode, expected="either 'same', 'full', or 'valid'")

    files_contain_data_cubes = fits.getdata(params.inFiles[0]).ndim == 3
    pad_vectors, reference_image_pad_vector = get_pad_vectors(shifts=shifts,
                                                            cube_mode=files_contain_data_cubes,
                                                            return_reference_image_pad_vector=True)

    # Assert that there are the same number of inFiles and psfFiles, which
    # should be the case after running the holography function.
    if not len(params.inFiles) == len(params.psfFiles):
        raise RuntimeError("The number of input files ({}) and PSF files ({}) do not match!".format(len(params.inFiles), len(params.psfFiles)))

    # Padding and Fourier transforming the images
    logger.info("Padding the images and PSFs...")
    for file_index, image_file in enumerate(params.inFiles):
        # Initialization
        image_pad_vector = pad_vectors[file_index]
        image_pad_vector.pop(0)
        if file_index == 0:
            # Get final image size
            img = fits.getdata(image_file)[0] # Remove time axis padding
            print("\tPadding data from {}".format(image_file))
            img = pad_array(array=img,
                            pad_vector=image_pad_vector,
                            mode=mode,
                            reference_image_pad_vector=reference_image_pad_vector)
            print('\tShift:', shifts[file_index])
            print('\tShape:', img.shape)

            # Get pad vector for PSFs
            psf_file = params.psfFiles[file_index]
            psf = fits.getdata(psf_file)[0]
            # Pad the Fpsf cube to have the same xz-extent as Fimg
            print("\tPadding data from {}".format(psf_file))
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
                raise ValueError("The Fourier transformed images and psfs have different shape, {} and {}. Something went wrong with the padding!".format(img.shape, psf.shape))

            # Initialize the enumerator and denominator
            enumerator = np.zeros(img.shape, dtype='complex128')
            denominator = np.zeros(img.shape, dtype='complex128')

        # Open PSF file
        print("\r\tFourier transforming image and PSF file {:4}/{:4}".format(file_index + 1, len(params.inFiles)), end='')
        psf_cube = fits.getdata(params.psfFiles[file_index])
        for frame_index, frame in enumerate(fits.getdata(image_file)):
            # Padding and transforming the image
            img = pad_array(array=frame,
                            pad_vector=pad_vectors[file_index],
                            mode=mode,
                            reference_image_pad_vector=reference_image_pad_vector)
            Fimg = fftshift(fft2(img))

            # Padding and Fourier transforming PSF
            psf = psf_cube[frame_index]
            psf = np.pad(psf, psf_pad_vector, mode='constant',)
            Fpsf = fftshift(fft2(psf))

            # Adding for the average
            enumerator += np.multiply(Fimg, np.conjugate(Fpsf))
            denominator += np.abs(np.square(Fpsf))

    print()

    # Compute the object.
    # Note that by this division implicitly does averaging. By this implicit
    # summing up of enumerator and denominator, this computation is cheaper
    # in terms of memory usage
    Fobject = np.divide(enumerator, denominator)

    return Fobject
