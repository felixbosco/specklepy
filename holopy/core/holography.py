import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile
from holopy.core.alignment import compute_shifts
from holopy.core.aperture import Aperture
from holopy.core.apodization import apodize
from holopy.core.ssa import ssa
# from holopy.algorithms.matchstarlist import MatchStarList
from holopy.algorithms.psfextraction import PSFExtraction
from holopy.algorithms.sourceextraction import SourceExtraction
# from holopy.algorithms.ssa import SSAReconstruction
from holopy.utils.plot import imshow
from holopy.utils.transferfunctions import otf



def holography(params, show=False, debug=False):
    """Execute the holographic image reconstruction following the algorithm
    outlined in Schoedel et al (2013, Section 3).
    """
    logging.info("Starting holographic reconstruction of {} files...".format(len(params.inFiles)))

    # (i-ii) Align cubes
    shifts = compute_shifts(files=params.inFiles,
                                    reference_file=params.alignmentReferenceFile,
                                    reference_file_index=0,
                                    lazy_mode=True,
                                    return_image_shape=False,
                                    debug=debug)
    pad_vectors = len(params.inFiles) * [((0, 0), (0, 0), (0, 0))]


    # (iii) Compute SSA reconstruction
    image = ssa(params.inFiles, outfile=params.outFile, tmp_dir=params.tmpDir)
    total_flux = np.sum(image) # Stored for flux conservation


    # Start iteration from steps (iv) thorugh (xi)
    while True:
        # (iv) Astrometry and photometry, i.e. StarFinder
        finder = SourceExtraction()
        finder.find_sources(image=image, starfinder_fwhm=params.starfinderFwhm, noise_threshold=params.noiseThreshold,
            background_subtraction=True, verbose=False)
        finder.writeto(params.allStarsFile)

        # (v) Select reference stars
        print("\tPlease copy your desired reference stars from the all stars file into the reference star file!")
        input("\tWhen you are done, just hit a key.")

        # (vi) PSF extraction
        algorithm = PSFExtraction(params)
        algorithm.extract(file_shifts=shifts, inspect_aperture=False)
        logging.info("Saved the extracted PSFs...")

        # (vii) Noise thresholding
        for file in params.psfFiles:
            with fits.open(file, mode='update') as hdulist:
                numberFrames = hdulist[0].header['NAXIS3']
                if not hasattr(params, 'psfNoiseMask'):
                    params.psfNoiseMask = _generate_noise_mask(hdulist[0].data[0], noise_reference_margin=params.noiseReferenceMargin)
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
        pass

        # (ix) Estimate object, following Eq. 1 (Schoedel et al., 2013)
        Fobject = evaluate_object(params, pad_vectors=pad_vectors)

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
        answer = input("\tDo you want to continue with one more iteration? [yes/no]")
        if answer.lower() in ['n', 'no']:
            break

    # Finally return the image
    return image





def _generate_noise_mask(frame, noise_reference_margin):
    """Estimate background (np.mean) and noise level (np.std) for every psf
    in a given outer annulus."""
    center = int((frame.shape[0] - 1) / 2)
    radius = center - noise_reference_margin
    tmp = Aperture(center, center, radius, frame, subset_only=False)
    return np.logical_not(tmp.data.mask)




def evaluate_object(params, pad_vectors, mode='same'):
    logging.info("Fourier transforming the images...")
    for index, file in enumerate(params.inFiles):
        img = fits.getdata(file)
        # logging.warning("Images are not expanded to the fill field of view.")
        print(pad_vectors)
        pad_vector = pad_vectors[index]
        img = np.pad(img, pad_vector)
        if mode == 'same':
            # Only in same mode, remove the margin outside the field of view
            # of the reference image, see align_cubes() method.
            img = img[:, pad_vector[1][0] : - pad_vector[1][1] -1, pad_vector[2][0]: - pad_vector[2][1] - 1]

        if index == 0:
            Fimg = fftshift(fft2(img))
        else:
            Fimg = np.concatenate((Fimg, fftshift(fft2(img))))


    # Clear memory
    img_shape = img.shape
    del img

    logging.info("Fourier transforming the PSFs...")
    for index, file in enumerate(params.psfFiles):
        psf = fits.getdata(file)
        if index == 0:
            # Pad the Fpsf cube to have the same xz-extent as Fimg
            print("\tPadding the PSFs")
            print('\tImages', img_shape)
            print('\tPSFs', psf.shape)
            dx = img_shape[1] - psf.shape[1]
            dy = img_shape[2] - psf.shape[2]

            pad_vector = ((0, 0), (int(np.floor(dx/2)), int(np.ceil(dx/2))), (int(np.floor(dy/2)), int(np.ceil(dy/2))))
            print('\tPad_width', pad_vector)
            psf = np.pad(psf, pad_vector)
            try:
                assert img_shape == psf.shape
            except:
                raise ValueError("The Fourier transformed images and psfs have different shape, {} and {}. Something went wrong with the padding!".format(Fimg_shape, Fpsf.shape))

            # Initialize Fpsf by the transforming the first cube
            Fpsf = fftshift(fft2(psf))
        else:
            Fpsf = np.concatenate((Fpsf, fftshift(fft2(np.pad(psf, pad_vector)))))

    # Clear memory
    del psf

    # Compute the object
    enumerator = np.mean(np.multiply(Fimg, np.conjugate(Fpsf)), axis=0)
    denominator = np.mean(np.abs(np.square(Fpsf)), axis=0)
    # denominator = np.ma.masked_values(denominator, 0.0)
    Fobject  = np.divide(enumerator, denominator)

    return Fobject
