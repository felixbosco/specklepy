import numpy as np

from astropy.io import fits

from specklepy.core.aperture import Aperture
from specklepy.core.fourierobject import FourierObject
from specklepy.core.psfextraction import ReferenceStars
from specklepy.core.reconstruction import Reconstruction
from specklepy.core.sourceextraction import extract_sources
from specklepy.io.filearchive import FileArchive
from specklepy.io.reconstructionfile import ReconstructionFile
from specklepy.exceptions import SpecklepyValueError
from specklepy.logging import logger
from specklepy.plotting.utils import imshow


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
    file_archive = FileArchive(file_list=params['PATHS']['inDir'], cards=[], dtypes=[])
    in_files = file_archive.files
    in_dir = file_archive.in_dir
    tmp_dir = params['PATHS']['tmpDir']

    # Input check
    if mode not in ['same', 'full', 'valid']:
        raise SpecklepyValueError('holography()', argname='mode', argvalue=mode,
                                  expected="either 'same', 'full', or 'valid'")

    if 'apodizationType' in params['APODIZATION']:
        # Catch deprecated parameter name
        logger.warning("Parameter 'apodizationType' is deprecated. Use 'type' instead!")
        params['APODIZATION']['type'] = params['APODIZATION']['apodizationType']
    if 'apodizationWidth' in params['APODIZATION']:
        # Catch deprecated parameter name
        logger.warning("Parameter 'apodizationWidth' is deprecated. Use 'radius' instead!")
        params['APODIZATION']['radius'] = params['APODIZATION']['apodizationWidth']
    if params['APODIZATION']['type'] is None or params['APODIZATION']['type'].lower() not in ['gaussian', 'airy']:
        logger.error(f"Apodization type has not been set or of wrong type ({params['APODIZATION']['type']})")
    if params['APODIZATION']['radius'] is None or not isinstance(params['APODIZATION']['radius'], (int, float)):
        logger.error(f"Apodization radius has not been set or of wrong type ({params['APODIZATION']['radius']})")

    # Initialize the outfile
    out_file = ReconstructionFile(filename=params['PATHS']['outFile'], files=in_files,
                                  cards={"RECONSTRUCTION": "Holography"}, in_dir=in_dir)

    # Initialize reconstruction
    reconstruction = Reconstruction(in_files=in_files, mode=mode, integration_method='ssa',
                                    reference_file=params['PATHS']['alignmentReferenceFile'],
                                    in_dir=in_dir, tmp_dir=tmp_dir, out_file=params['PATHS']['outFile'],
                                    var_ext=params['OPTIONS']['varianceExtensionName'],
                                    box_indexes=params['OPTIONS']['box_indexes'], debug=debug)

    # (i-ii) Align cubes
    reconstruction.align_cubes()
    shifts = reconstruction.shifts

    # (iii) Compute SSA reconstruction
    image = reconstruction.coadd_long_exposures(save=True)
    if isinstance(image, tuple):
        # SSA returned a reconstruction image and a variance image
        image, image_var = image
    total_flux = np.sum(image)  # Stored for flux conservation

    # Start iteration from steps (iv) through (xi)
    while True:
        # (iv) Astrometry and photometry, i.e. StarFinder
        # (v) Select reference stars
        extract_sources(image=image,
                        fwhm=params['STARFINDER']['starfinderFwhm'],
                        noise_threshold=params['STARFINDER']['noiseThreshold'],
                        background_subtraction=True,
                        write_to=params['PATHS']['allStarsFile'],
                        star_finder='DAO', select=params['PATHS']['refSourceFile'], debug=debug)


        # print("\tPlease copy your desired reference stars from the all stars file into the reference star file!")
        # input("\tWhen you are done, hit a ENTER.")

        # (vi) PSF extraction
        psf_stars = ReferenceStars(psf_radius=params['PSFEXTRACTION']['psfRadius'],
                                   reference_source_file=params['PATHS']['refSourceFile'], in_files=in_files,
                                   save_dir=tmp_dir, in_dir=in_dir,
                                   field_segmentation=params['PSFEXTRACTION']['fieldSegmentation'])
        if params['PSFEXTRACTION']['mode'].lower() == 'epsf':
            psf_files = psf_stars.extract_epsfs(file_shifts=shifts, debug=debug)
        elif params['PSFEXTRACTION']['mode'].lower() in ['mean', 'median', 'weighted_mean']:
            psf_files = psf_stars.extract_psfs(file_shifts=shifts, mode=params['PSFEXTRACTION']['mode'].lower(),
                                               debug=debug)
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
        f_object.apodize(type=params['APODIZATION']['type'], radius=params['APODIZATION']['radius'])

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
    extract_sources(image=image, fwhm=params['STARFINDER']['starfinderFwhm'],
                    noise_threshold=params['STARFINDER']['noiseThreshold'], background_subtraction=True,
                    write_to=params['PATHS']['allStarsFile'], star_finder='DAO', debug=debug)

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
