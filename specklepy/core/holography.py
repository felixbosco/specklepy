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


def holography(params, debug=False):
    """Execute the holographic image reconstruction.

    The holographic image reconstruction is an algorithm as outlined, eg. by Schoedel et al (2013, Section 3). This
    function follows that algorithm, see comments in the code. Most of the important functions are imported from other
    modules of specklepy.

    Args:
        params (dict):
            Dictionary that carries all important parameters.
        debug (bool, optional):
            Set to True to inspect intermediate results.
            Default is False.

    Returns:
        image (np.ndarray): The image reconstruction.
    """

    logger.info(f"Starting holographic reconstruction...")

    # Extract individual parameter dictionaries
    paths = params.get('PATHS', {})
    alignment = params.get('ALIGNMENT', {})
    apodization = params.get('APODIZATION', {})
    psf_extraction = params.get('PSFEXTRACTION', {})
    source_extraction = params.get('STARFINDER', {})

    # Create file archive and export paths to be directly available
    file_archive = FileArchive(file_list=paths.get('inDir'), cards=[], dtypes=[],
                               table_format=paths.get('tableFormat', 'ascii.no_header'))
    in_files = file_archive.files
    in_dir = file_archive.in_dir
    tmp_dir = paths.get('tmpDir')

    # Check input mode
    mode = params['ALIGNMENT'].get('reconstructionMode')
    if mode not in ['same', 'full', 'valid']:
        raise SpecklepyValueError('holography()', argname='mode', argvalue=mode,
                                  expected="either 'same', 'full', or 'valid'")

    # Check apodization parameters
    if 'apodizationType' in apodization:
        # Catch deprecated parameter name
        logger.warning("Parameter 'apodizationType' is deprecated. Use 'type' instead!")
        apodization['type'] = apodization['apodizationType']
    if 'apodizationWidth' in apodization:
        # Catch deprecated parameter name
        logger.warning("Parameter 'apodizationWidth' is deprecated. Use 'radius' instead!")
        apodization['radius'] = apodization['apodizationWidth']
    if apodization['type'] is None or apodization['type'].lower() not in ['gaussian', 'airy']:
        logger.error(f"Apodization type has not been set or of wrong type ({apodization['type']})")
    if apodization['radius'] is None or not isinstance(apodization['radius'], (int, float)):
        logger.error(f"Apodization radius has not been set or of wrong type ({apodization['radius']})")

    # Initialize the outfile
    out_file = ReconstructionFile(filename=paths.get('outFile'), files=in_files,
                                  cards={"RECONSTRUCTION": "Holography"}, in_dir=in_dir)

    # Initialize reconstruction
    reconstruction = Reconstruction(in_files=in_files, mode=mode, integration_method='ssa',
                                    reference_file=alignment.get('referenceFile'),
                                    in_dir=in_dir, tmp_dir=tmp_dir, out_file=paths.get('outFile'),
                                    variance_extension=params['EXTNAMES']['varianceExtension'],
                                    box_indexes=None, custom_mask=paths.get('maskFile'), debug=debug)
    reconstruction.assert_dirs()

    # (i-ii) Align cubes
    # Compute the first alignment based on collapsed images (and variance images)
    mask_hot_pixels = alignment.get('maskHotPixels', False)
    reconstruction.align_cubes(integration_method='collapse', alignment_mode=alignment.get('mode', 'correlation'),
                               mask_hot_pixels=mask_hot_pixels)

    # Repeat alignment in SSA mode, if not requested otherwise
    if alignment.get('integrationMode', 'ssa') != 'collapse':
        if psf_extraction.get('psfRadius') is not None:
            reconstruction.select_box(radius=psf_extraction.get('psfRadius'))
        reconstruction.long_exp_files = reconstruction.create_long_exposures(integration_method='ssa',
                                                                             mask_hot_pixels=mask_hot_pixels,
                                                                             shifts=reconstruction.shifts)

        reconstruction.align_cubes(integration_method='ssa', alignment_mode=alignment.get('mode', 'correlation'),
                                   mask_hot_pixels=mask_hot_pixels)
    shifts = reconstruction.shifts

    # (iii) Compute SSA reconstruction
    image, image_var = reconstruction.coadd_long_exposures(save=True)
    total_flux = np.sum(image)  # Stored for flux conservation

    # Save SSA reconstruction to serparate file, if requested
    if paths.get('ssaFile') is not None:
        ssa_file = ReconstructionFile(filename=paths.get('ssaFile'), files=in_files, in_dir=in_dir)
        ssa_file.data = image
        ssa_file.new_extension(name=params['EXTNAMES']['varianceExtension'], data=image_var)

    # Start iteration from steps (iv) through (xi)
    while True:
        # (iv) Astrometry and photometry, i.e. StarFinder
        # (v) Select reference stars
        if psf_extraction.get('select', 'gui') == 'gui':
            extract_sources(image=image,
                            fwhm=source_extraction.get('starfinderFwhm'),
                            noise_threshold=source_extraction.get('noiseThreshold'),
                            background_subtraction=True,
                            write_to=paths.get('allStarsFile'),
                            star_finder=source_extraction.get('algorithm'),
                            select=paths.get('refSourceFile'), debug=debug)
        else:
            logger.info("\tPlease copy the desired reference stars from the all stars file into the reference star "
                        "file!")
            input("\tWhen you are done, hit ENTER.")

        # (vi) PSF extraction
        psf_stars = ReferenceStars(psf_radius=psf_extraction.get('psfRadius'),
                                   reference_source_file=paths.get('refSourceFile'), in_files=in_files,
                                   save_dir=tmp_dir, in_dir=in_dir,
                                   field_segmentation=psf_extraction.get('fieldSegmentation'))
        if psf_extraction.get('mode').lower() == 'epsf':
            psf_files = psf_stars.extract_epsfs(file_shifts=shifts, debug=debug)
        elif psf_extraction.get('mode').lower() in ['mean', 'median', 'weighted_mean']:
            psf_files = psf_stars.extract_psfs(file_shifts=shifts, mode=psf_extraction.get('mode').lower(),
                                               debug=debug)
        else:
            raise RuntimeError(f"PSF extraction mode '{psf_extraction.get('mode')}' is not understood!")
        logger.info("Saved the extracted PSFs...")

        # (vii) Noise thresholding
        psf_noise_mask = None
        for file in psf_files:
            with fits.open(file, mode='update') as hdu_list:
                n_frames = hdu_list[0].header['NAXIS3']
                if psf_noise_mask is None:
                    psf_noise_mask = get_noise_mask(hdu_list[0].data[0],
                                                    noise_reference_margin=
                                                    psf_extraction.get('noiseReferenceMargin'))
                for index in range(n_frames):
                    reference = np.ma.masked_array(hdu_list[0].data[index], mask=psf_noise_mask)
                    background = np.mean(reference)
                    noise = np.std(reference)
                    update = np.maximum(hdu_list[0].data[index] - background -
                                        psf_extraction.get('noiseThreshold') * noise, 0.0)
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
        f_object = FourierObject(in_files, psf_files, shifts=shifts, mode=mode, in_dir=in_dir,
                                 custom_mask=paths.get('maskFile'))
        f_object.coadd_fft(mask_hot_pixels=alignment.get('maskHotPixels', False),
                           bootstrap=params.get('BOOTSTRAP').get('numberImages'))

        # (x) Apodization
        f_object.apodize(type=apodization.get('type'), radius=apodization.get('radius'))

        # (xi) Inverse Fourier transform to retain the reconstructed image
        image, bootstrap_images = f_object.ifft(total_flux=total_flux)

        # Inspect the latest reconstruction
        if debug:
            imshow(image)

        # Save the latest reconstruction image to outfile
        out_file.data = image

        # Compute uncertainty from bootstrap reconstructions
        if bootstrap_images is not None:
            var = np.var(np.array(bootstrap_images), axis=0)
            out_file.update_extension(params['EXTNAMES']['varianceExtension'], var)
            if params['BOOTSTRAP'].get('saveImages', False):
                logger.info("Saving bootstrap images...")
                bs_file = ReconstructionFile(filename=paths.get('outFile').replace('.fits', '_bs.fits'), files=in_files,
                                             cards={"RECONSTRUCTION": "Holography"}, in_dir=in_dir,
                                             shape=(len(bootstrap_images), image.shape[1], image.shape[2]))
                for b, bootstrap_image in enumerate(bootstrap_images):
                    bs_file.update_frame(frame_index=b, data=bootstrap_image)

        # Ask the user whether the iteration shall be continued or not
        answer = input("\tDo you want to continue with one more iteration? [yes/no]\n\t")
        if answer.lower() in ['n', 'no']:
            break

    # Repeat astrometry and photometry, i.e. StarFinder on final image
    fwhm = apodization.get('radius')
    if 'gauss' in apodization.get('type').lower():
        fwhm *= 2.35
    extract_sources(image=image, fwhm=fwhm, noise_threshold=source_extraction.get('noiseThreshold'),
                    background_subtraction=True, write_to=paths.get('allStarsFile'),
                    star_finder=source_extraction.get('algorithm'), debug=debug)

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
