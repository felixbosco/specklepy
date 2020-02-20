import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from specklepy.core import alignment
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.outfile import Outfile
from specklepy.io.reconstructionfile import ReconstructionFile
from specklepy.logging import logging



def ssa(files, mode='same', reference_file=None, outfile=None, tmp_dir=None, lazy_mode=True, debug=False, **kwargs):
    """Compute the SSA reconstruction of a list of files.

    The simple shift-and-add (SSA) algorithm makes use of the structure of
    typical speckle patterns, i.e. short-exposure point-spread functions (PSFs).
    These show multiple peaks resembling the diffraction-limited PSF of coherent
    fractions within the telescope aperture. Under good conditions or on small
    telescopes, there is typically one largest coherent atmospheric cell and
    therefore, speckle PSFs typically show one major intensity peak. The
    algorithm makes use of this fact and identifies the emission peak in a given
    observation frame, assuming that this always belongs to the same star, and
    aligns all frames on the coordinate of the emission peak.

    See Bates & Cady (1980) for references.

    Args:
        files (list):
            List of complete paths to the fits files that shall be considered
            for the SSA reconstruction.
        mode (str):
            Name of the reconstruction mode: In 'same' mode, the reconstruction
            covers the same field of view of the reference file. In 'full' mode,
            every patch of the sky that is covered by at least one frame will be
            contained in the final reconstruction.
        reference_file (str, int, optional):
            Path to a reference file or index of the file in files, relative to
            which the shifts are computed. See
            specklepy.core.aligment.get_shifts for details. Default is 0.
        outfile (specklepy.io.recfile, optional): Object to write the result to,
            if provided.
        tmp_dir (str, optional):
            Path of a directory in which the temporary results are stored in.
        lazy_mode (bool, optional):
            Set to False, to enforce the alignment of a single file with respect
            to the reference file. Default is True.
        debug (bool, optional): Set to True to inspect intermediate results.
            Default is False.

    Returns:
        reconstruction (np.ndarray):
            The image reconstruction. The size depends on the mode argument.
    """

    logging.info("Starting SSA reconstruction...")
    # Check parameters
    if not isinstance(files, list):
        if isinstance(files, str):
            files = [files]
        else:
            raise SpecklepyTypeError('ssa()', argname='files', argtype=type(files), expected='list')

    if isinstance(mode, str):
        if not mode in ['same', 'full', 'valid']:
            raise SpecklepyValueError('ssa()', argname='mode', argvalue=mode, expected="'same', 'full' or 'valid'")
    else:
        raise SpecklepyTypeError('ssa()', argname='mode', argtype=type(mode), expected='str')

    if reference_file is None:
        reference_file = files[0]
    elif isinstance(reference_file, int):
        reference_file = files[reference_file]
    elif not isinstance(reference_file, str):
        raise SpecklepyTypeError('ssa()', argname='reference_file', argtype=type(reference_file), expected='str or int')

    if outfile is None:
        pass
    elif isinstance(outfile, str):
        outfile = ReconstructionFile(files=files, filename=outfile, cards={"RECONSTRUCTION": "SSA"})
    elif isinstance(outfile, ReconstructionFile):
        pass
    else:
        raise SpecklepyTypeError('ssa()', argname='outfile', argtype=type(outfile), expected='str')

    if tmp_dir is not None:
        if isinstance(tmp_dir, str) and not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)

    if not isinstance(lazy_mode, bool):
        raise SpecklepyTypeError('ssa()', argname='lazy_mode', argtype=type(lazy_mode), expected='bool')

    if 'variance_extension_name' in kwargs.keys():
        var_ext = kwargs['variance_extension_name']
    else:
        var_ext = 'VAR'

    # Do not align just a single file
    if lazy_mode and len(files) == 1:
        with fits.open(files[0]) as hdulist:
            cube = hdulist[0].data
            try:
                var_cube = hdulist[var_ext].data
                reconstruction, reconstruction_var = coadd_frames(cube, var_cube=var_cube)
            except:
                reconstruction = coadd_frames(cube)

    # Align reconstructions if multiple files are given
    else:
        # Compute temporary reconstructions of the individual cubes
        tmp_files = []
        for index, file in enumerate(files):
            with fits.open(files[0]) as hdulist:
                cube = hdulist[0].data
                try:
                    var_cube = hdulist[var_ext].data
                    tmp, tmp_var = coadd_frames(cube, var_cube=var_cube)
                except:
                    tmp = coadd_frames(cube)
            tmp_file = os.path.basename(file).replace(".fits", "_ssa.fits")
            tmp_file = os.path.join(tmp_dir, tmp_file)
            logging.info("Saving interim SSA reconstruction of cube to {}".format(tmp_file))
            tmp_file_object = Outfile(tmp_file, shape=tmp.shape, verbose=False)
            tmp_file_object.data = tmp
            if 'tmp_var' in locals():
                tmp_file_object.new_extension(var_ext, data=tmp_var)
            # fits.writeto(tmp_file, tmp, overwrite=True)
            tmp_files.append(tmp_file)

        # Align tmp reconstructions and add up
        file_shifts, image_shape = alignment.get_shifts(tmp_files, reference_file=reference_file, return_image_shape=True, lazy_mode=True)
        pad_vectors, ref_pad_vector = alignment.get_pad_vectors(file_shifts, cube_mode=(len(image_shape) == 3), return_reference_image_pad_vector=True)
        # reconstruction = np.zeros(image_shape)
        for index, file in enumerate(tmp_files):
            with fits.open(file) as hdulist:
                tmp_image = hdulist[0].data
                try:
                    tmp_image_var = hdulist[var_ext].data
                except:
                    # No VAR extension in file
                    pass
            if index is 0:
                reconstruction = alignment.pad_array(tmp_image, pad_vectors[index], mode=mode, reference_image_pad_vector=ref_pad_vector)
                try:
                    reconstruction_var = alignment.pad_array(tmp_image_var, pad_vectors[index], mode=mode, reference_image_pad_vector=ref_pad_vector)
                except:
                    pass
            else:
                reconstruction += alignment.pad_array(tmp_image, pad_vectors[index], mode=mode, reference_image_pad_vector=ref_pad_vector)
                try:
                    reconstruction_var += alignment.pad_array(tmp_image_var, pad_vectors[index], mode=mode, reference_image_pad_vector=ref_pad_vector)
                except:
                    pass

    logging.info("Reconstruction finished...")

    # Save the result to an Outfile
    if outfile is not None:
        outfile.data = reconstruction
        if 'reconstruction_var' in locals():
            outfile.new_extension(name=var_ext, data=reconstruction_var)

    if 'reconstruction_var' in locals():
        return reconstruction, reconstruction_var
        # Otherwise:
    return reconstruction



def coadd_frames(cube, var_cube=None):
    """Compute the simple shift-and-add (SSA) reconstruction of a data cube.

    This function uses the SSA algorithm to coadd frames of a cube. If
    provided, this function coadds the variances within a var cube considering
    the exact same shifts.

    Args:
        cube (np.ndarray, ndim=3):
            Data cube which is integrated along the zero-th axis.
        var_cube (np.ndarray, ndim=3, optional):
            Data cube of variances which is integrated along the zero-th axis
            with the same shifts as the cube.
    Returns:
        coadded (np.ndarray, ndim=2):
            SSA-integrated frames of the input cube.
    """

    if not isinstance(cube, np.ndarray):
        raise SpecklepyTypeError('coadd_frames()', argname='cube', argtype=type(cube), expected='np.ndarray')
    if cube.ndim is not 3:
        raise SpecklepyValueError('coadd_frames()', argname='cube.ndim', argvalue=cube.ndim, expected='3')

    if var_cube is not None and not isinstance(var_cube, np.ndarray):
        raise SpecklepyTypeError('coadd_frames()', argname='var_cube', argtype=type(var_cube), expected='np.ndarray')
    if var_cube is not None and var_cube.shape != cube.shape:
        raise SpecklepyValueError('coadd_frames()', argname='var_cube.shape', argvalue=var_cube.shape, expected=str(cube.shape))

    # Compute shifts
    peak_indizes = np.zeros((cube.shape[0], 2), dtype=int)
    for index, frame in enumerate(cube):
        peak_indizes[index] = np.array(np.unravel_index(np.argmax(frame, axis=None), frame.shape), dtype=int)


    # Compute shifts from indizes
    peak_indizes = peak_indizes.transpose()
    xmean, ymean = np.mean(np.array(peak_indizes), axis=1)
    xmean = int(xmean)
    ymean = int(ymean)
    shifts = np.array([xmean - peak_indizes[0], ymean - peak_indizes[1]])
    shifts = shifts.transpose()

    # Shift frames and add to coadded
    coadded = np.zeros(cube[0].shape)
    pad_vectors, ref_pad_vector = alignment.get_pad_vectors(shifts, cube_mode=False, return_reference_image_pad_vector=True)
    for index, frame in enumerate(cube):
        coadded += alignment.pad_array(frame, pad_vectors[index], mode='same', reference_image_pad_vector=ref_pad_vector)

    if var_cube is not None:
        var_coadded = np.zeros(coadded.shape)
        for index, frame in enumerate(var_cube):
            var_coadded += alignment.pad_array(frame, pad_vectors[index], mode='same', reference_image_pad_vector=ref_pad_vector)
        return coadded, var_coadded

    return coadded

