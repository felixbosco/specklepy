import os

from specklepy.core import Reconstruction
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io import FileArchive
from specklepy.logging import logger


def ssa(files, mode='same', reference_file=None, outfile=None, in_dir=None, tmp_dir=None, variance_extension=None,
        aperture_radius=None, integration_method='ssa', alignment_method='correlation', mask_hot_pixels=False,
        mask_file=None, debug=False, **kwargs):
    """Compute the SSA reconstruction of a list of files.

    The simple shift-and-add (SSA) algorithm makes use of the structure of typical speckle patterns, i.e.
    short-exposure point-spread functions (PSFs). These show multiple peaks resembling the diffraction-limited PSF of
    coherent fractions within the telescope aperture. Under good conditions or on small telescopes, there is typically
    one largest coherent atmospheric cell and therefore, speckle PSFs typically show one major intensity peak. The
    algorithm makes use of this fact and identifies the emission peak in a given observation frame, assuming that this
    always belongs to the same star, and aligns all frames on the coordinate of the emission peak.

    See Bates & Cady (1980) for references.

    Args:
        files (list or array_like):
            List of complete paths to the fits files that shall be considered for the SSA reconstruction.
        mode (str):
            Name of the reconstruction mode: In 'same' mode, the reconstruction covers the same field of view of the
            reference file. In 'full' mode, every patch of the sky that is covered by at least one frame will be
            contained in the final reconstruction.
        reference_file (str, int, optional):
            Path to a reference file or index of the file in files, relative to which the shifts are computed. See
            specklepy.core.aligment.get_shifts for details. Default is 0.
        outfile (str, optional):
            Object to write the result to, if provided.
        in_dir (str, optional):
            Path to the files. `None` is substituted by an empty string.
        tmp_dir (str, optional):
            Path of a directory in which the temporary results are stored in.
        variance_extension (str, optional):
            Name of the variance extension in the FITS files. Falls back to `'VAR'`.
        aperture_radius (list, optional):
            Constraining the search for the intensity peak to an aperture with the specified radius. The position of
            the aperture is selected graphically. The code is estimating the intensity peak within the full frames if
            not provided.
        integration_method (str, optional):
            Method for creating individual long exposures. Using SSA by default, but can also use 'collapse' for a
            straight integration along the time axis (in case of faint reference sources).
        alignment_method (str, optional):
            Method for aligning cubes: Can be either 'correlation', 'peak', or 'sources'
        mask_hot_pixels (bool, optional):
            Mask hot pixels prior to alignment.
        mask_file (str, optional):
            Name of a file containing a mask, which will be applied additionally to the files own mask extension.
        debug (bool, optional):
            Show debugging information. Default is False.

    Returns:
        reconstruction (np.ndarray):
            The image reconstruction. The size depends on the mode argument.
        reconstruction_var (np.ndarray):
            The variance map of the reconstructed image. The size is the same as `reconstruction`.
    """

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')
        logger.handlers[0].setLevel('DEBUG')
        logger.info("Set logging level to DEBUG")

    # Check parameters
    file_archive = FileArchive(file_list=files, table_format=kwargs.get('tableFormat', 'ascii.no_header'))
    files = file_archive.files
    in_dir = file_archive.file_path

    if isinstance(mode, str):
        if mode not in ['same', 'full', 'valid']:
            raise SpecklepyValueError('ssa()', argname='mode', argvalue=mode, expected="'same', 'full' or 'valid'")
    else:
        raise SpecklepyTypeError('ssa()', argname='mode', argtype=type(mode), expected='str')

    if outfile is None or isinstance(outfile, str):
        pass
    else:
        raise SpecklepyTypeError('ssa()', argname='outfile', argtype=type(outfile), expected='str')

    if tmp_dir is not None:
        if isinstance(tmp_dir, str) and not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)

    if variance_extension is None:
        variance_extension = 'VAR'

    # Initialize the reconstruction
    logger.info("Starting SSA reconstruction...")
    reconstruction = Reconstruction(in_files=files, mode=mode, integration_method=integration_method,
                                    reference_file=reference_file,
                                    in_dir=in_dir, tmp_dir=tmp_dir, out_file=outfile,
                                    variance_extension=variance_extension,
                                    box_indexes=None, custom_mask=mask_file, debug=debug)

    # Compute the first alignment based on collapsed images (and variance images)
    reconstruction.align_cubes(integration_method='collapse', alignment_mode=alignment_method,
                               mask_hot_pixels=mask_hot_pixels)

    # Repeat alignment in SSA mode, if not requested otherwise
    if integration_method != 'collapse':
        if aperture_radius is not None:
            reconstruction.select_box(radius=aperture_radius)
        reconstruction.long_exp_files = reconstruction.create_long_exposures(integration_method='ssa',
                                                                             mask_hot_pixels=mask_hot_pixels,
                                                                             shifts=reconstruction.alignment.shifts)

        reconstruction.align_cubes(integration_method=integration_method, alignment_mode=alignment_method,
                                   mask_hot_pixels=mask_hot_pixels)
    reconstruction_image, reconstruction_var = reconstruction.coadd_long_exposures(save=True)

    # Return reconstruction and the variance map, which may be `None` type
    return reconstruction_image, reconstruction_var
