import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
import os
import glob
from datetime import datetime
from astropy.io import fits

from specklepy.logging import logging
from specklepy.io.parameterset import ParameterSet
from specklepy.io.outfile import Outfile
from specklepy.io.recfile import RECfile
from specklepy.core.alignment import get_shifts, get_pad_vectors, pad_array
from specklepy.core.aperture import Aperture
from specklepy.core.apodization import apodize
from specklepy.core.ssa import ssa
from specklepy.core.psfextraction import ReferenceStars
from specklepy.core.sourceextraction import find_sources
from specklepy.utils.plot import imshow
from specklepy.utils.transferfunctions import otf



def holography(params, mode='same', debug=False):
    """Execute the holographic image reconstruction following the algorithm
    outlined in Schoedel et al (2013, Section 3).

    Long description ...

    Args:
        params (specklepy.io.parameterset.ParameterSet):
        mode (str, optional): Define the size of the output image as 'same'
            to the reference image or expanding to include the 'full'
            covered field. Default is 'same'.
        debug (bool, optional): Set to True to inspect intermediate results.
            Default is False.

    Returns:
        image (np.ndarray): The image reconstruction.
    """

    logging.info("Starting holographic reconstruction of {} files...".format(len(params.inFiles)))
    if not isinstance(params, ParameterSet):
        logging.warn("specklepy.core.holography.holography received params argument \
                        of type <{}> instead of the expected type \
                        specklepy.io.parameterset.ParameterSet. This may cause \
                        unforeseen errors.".format(type(params)))
    if mode not in ['same', 'full', 'valid']:
        raise ValueError("specklepy.core.holography.holography received mode \
                            argument '{}', but must be either 'same', 'full', \
                            or 'valid'.".format(mode))

    # Initialize the outfile
    params.outFile = RECfile(filename=params.outFile, files=params.inFiles, cards={"RECONSTRUCTION": "Holography"})

    # (i-ii) Align cubes
    shifts = get_shifts(files=params.inFiles,
                            reference_file=params.alignmentReferenceFile,
                            lazy_mode=True,
                            return_image_shape=False,
                            debug=debug)

    # (iii) Compute SSA reconstruction
    image = ssa(params.inFiles, outfile=params.outFile, tmp_dir=params.tmpDir)
    total_flux = np.sum(image) # Stored for flux conservation

    # Start iteration from steps (iv) thorugh (xi)
    while True:
        # (iv) Astrometry and photometry, i.e. StarFinder
        find_sources(image=image, fwhm=params.starfinderFwhm, noise_threshold=params.signalToNoiseThreshold,
            background_subtraction=False, writeto=params.allStarsFile, starfinder='DAO', verbose=False)

        # (v) Select reference stars
        print("\tPlease copy your desired reference stars from the all stars file into the reference star file!")
        input("\tWhen you are done, just hit a key.")

        # (vi) PSF extraction
        refStars = ReferenceStars(params)
        refStars.extract_psfs(file_shifts=shifts, debug=False)
        logging.info("Saved the extracted PSFs...")

        # (vii) Noise thresholding
        for file in params.psfFiles:
            with fits.open(file, mode='update') as hdulist:
                numberFrames = hdulist[0].header['NAXIS3']
                if not hasattr(params, 'psfNoiseMask'):
                    params.psfNoiseMask = get_noise_mask(hdulist[0].data[0], noise_reference_margin=params.noiseReferenceMargin)
                for index in range(numberFrames):
                    reference = np.ma.masked_array(hdulist[0].data[index], mask=params.psfNoiseMask)
                    background = np.mean(reference)
                    noise = np.std(reference)
                    update = np.maximum(hdulist[0].data[index] - background - params.noiseThreshold * noise, 0.0)
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
        logging.info("Apodizing the object...")
        Fobject = apodize(Fobject, params.apodizationType, radius=params.apodizationWidth)

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

    # Finally return the image
    return image





def get_noise_mask(frame, noise_reference_margin):
    """Create an annulus-like mask within a given aperture for measuring noise
    and (sky) background.

    Args:
        frame ():
        noise_reference_margin (int): Width of the reference annulus:

    Returns:
        annulus_mask (np.ndarray, dtype=bool):
    """
    center = int((frame.shape[0] - 1) / 2)
    radius = center - noise_reference_margin
    tmp = Aperture(center, center, radius, data=frame, crop=False)
    annulus_mask = np.logical_not(tmp.data.mask)
    return annulus_mask




def get_Fourier_object(params, shifts, mode='same'):
    """Reconstruction of the Fourier transformed object with Eq. 1 (Schoedel
    et al., 2013).

    Long description...

    Args:
        params (ParameterSet):
        shifts (list):
        mode (str, optional): Define the size of the output image as 'same'
            to the reference image or expanding to include the 'full'
            covered field. Default is 'same'.

    Returns:
        Fobject (np.ndarray, dtype=complex128): Fourier transformed object as a
            complex128 np.ndarray.
    """

    if not isinstance(params, ParameterSet):
        logging.warn("specklepy.core.holography.get_Fourier_object received params argument \
                        of type <{}> instead of the expected type \
                        specklepy.io.parameterset.ParameterSet. This may cause \
                        unforeseen errors.".format(type(params)))
    if mode not in ['same', 'full', 'valid']:
        raise ValueError("specklepy.core.holography.get_Fourier_object received mode \
                            argument '{}', but must be either 'same', 'full', \
                            or 'valid'.".format(mode))

    pad_vectors, reference_image_pad_vector = get_pad_vectors(shifts=shifts,
                                    array_shape=fits.getdata(params.inFiles[0]).shape,
                                    reference_image_shape=(1024, 1024),
                                    mode='same')

    # Assert that there are the same number of inFiles and psfFiles, which
    # should be the case after running the holography function.
    if not len(params.inFiles) == len(params.psfFiles):
        raise RuntimeError("The number of input files ({}) and PSF files ({}) do not match!".format(len(params.inFiles), len(params.psfFiles)))

    # Padding and Fourier transforming the images
    logging.info("Padding the images and PSFs...")
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
            psf = np.pad(psf, psf_pad_vector)
            try:
                assert img.shape == psf.shape
            except:
                raise ValueError("The Fourier transformed images and psfs have different shape, {} and {}. Something went wrong with the padding!".format(img.shape, psf.shape))

            # Initialize the enumerator and denominator
            enumerator = np.zeros(img.shape, dtype='complex128')
            denominator = np.zeros(img.shape, dtype='complex128')

        # Open PSF file
        print("\rFourier transforming image and PSF file {:4}/{:4}".format(file_index + 1, len(params.inFiles)), end='')
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
            psf = np.pad(psf, psf_pad_vector)
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
