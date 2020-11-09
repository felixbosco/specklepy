from IPython import embed
import numpy as np
import os

from specklepy.core.reconstruction import Reconstruction
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.logging import logger


def ssa(files, mode='same', reference_file=None, outfile=None, in_dir=None, tmp_dir=None, box_indexes=None,
        debug=False, **kwargs):
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
        outfile (specklepy.io.recfile, optional):
            Object to write the result to, if provided.
        in_dir (str, optional):
            Path to the files. `None` is substituted by an empty string.
        tmp_dir (str, optional):
            Path of a directory in which the temporary results are stored in.
        box_indexes (list, optional):
            Constraining the search for the intensity peak to the specified box. Searching the full frames if not
            provided.
        debug (bool, optional):
            Show debugging information. Default is False.

    Returns:
        reconstruction (np.ndarray):
            The image reconstruction. The size depends on the mode argument.
    """

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')
        logger.handlers[0].setLevel('DEBUG')
        logger.info("Set logging level to DEBUG")
    # Check parameters
    if not isinstance(files, (list, np.ndarray)):
        if isinstance(files, str):
            files = [files]
        else:
            raise SpecklepyTypeError('ssa()', argname='files', argtype=type(files), expected='list')

    if isinstance(mode, str):
        if mode not in ['same', 'full', 'valid']:
            raise SpecklepyValueError('ssa()', argname='mode', argvalue=mode, expected="'same', 'full' or 'valid'")
    else:
        raise SpecklepyTypeError('ssa()', argname='mode', argtype=type(mode), expected='str')

    if outfile is None or isinstance(outfile, str):
        pass
    else:
        raise SpecklepyTypeError('ssa()', argname='outfile', argtype=type(outfile), expected='str')

    if in_dir is None:
        in_dir = ''

    if tmp_dir is not None:
        if isinstance(tmp_dir, str) and not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)

    if 'variance_extension_name' in kwargs.keys():
        var_ext = kwargs['variance_extension_name']
    else:
        var_ext = 'VAR'

    # Initialize the reconstruction
    logger.info("Starting SSA reconstruction...")
    reconstruction = Reconstruction(in_files=files, mode=mode, integration_method='ssa',
                                    reference_file=reference_file,
                                    in_dir=in_dir, tmp_dir=tmp_dir, out_file=outfile,
                                    var_ext=var_ext,
                                    box_indexes=box_indexes, debug=debug)

    # Compute the aligned and co-added image (and variance image)
    reconstruction_image, reconstruction_var = reconstruction.coadd_long_exposures(save=True)

    # Return reconstruction (and the variance map if computed)
    if reconstruction_var is not None:
        return reconstruction_image, reconstruction_var
    return reconstruction_image
