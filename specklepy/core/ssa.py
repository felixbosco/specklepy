import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from specklepy.logging import logging
from specklepy.io.recfile import RECfile
from specklepy.core import alignment


def ssa(files, mode='same', reference_file=0, outfile=None, tmp_dir=None, lazy_mode=True, debug=False, **kwargs):
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
        files = [files]

    if isinstance(reference_file, int):
        reference_file = files[reference_file]
    elif not isinstance(reference_file, str):
        raise TypeError("The function get_shifts received reference_file argument of type {}, but needs be int or str, i.e. a file name.".format(type(reference_file)))

    if outfile is not None and not isinstance(outfile, RECfile):
        if isinstance(outfile, str):
            outfile = RECfile(files=files, filename=outfile, cards={"RECONSTRUCTION": "SSA"})
        else:
            raise TypeError("specklepy.core.ssa.ssa received outfile argument of wrong type <{}>!".format(type(outfile)))

    # Do not align just a single file
    if lazy_mode and len(files) == 1:
        cube = fits.getdata(files[0])
        reconstruction = coadd_frames(cube)

    # Align reconstructions if multiple files are given
    else:
        # Compute temporary reconstructions of the individual cubes
        tmp_files = []
        for index, file in enumerate(files):
            cube = fits.getdata(file)
            tmp = coadd_frames(cube)
            tmp_file = os.path.basename(file).replace(".fits", "_ssa.fits")
            tmp_file = os.path.join(tmp_dir, tmp_file)
            logging.info("Saving interim SSA reconstruction of cube to {}".format(tmp_file))
            fits.writeto(tmp_file, tmp, overwrite=True)
            tmp_files.append(tmp_file)

        # Align tmp reconstructions and add up
        file_shifts, image_shape = alignment.get_shifts(tmp_files, reference_file=reference_file, return_image_shape=True, lazy_mode=True)
        pad_vectors, ref_pad_vector = alignment.get_pad_vectors(file_shifts, image_shape, image_shape, mode='same')
        reconstruction = np.zeros(image_shape)
        for index, file in enumerate(tmp_files):
            tmp_image = fits.getdata(file)
            reconstruction += alignment.pad_array(tmp_image, pad_vectors[index], mode='same', reference_image_pad_vector=ref_pad_vector)

    logging.info("Reconstruction finished...")

    # Save the result to an Outfile
    if outfile is not None:
        outfile.data = reconstruction

    return reconstruction



def coadd_frames(cube):
    """
    Compute the simple shift-and-add (SSA) reconstruction of a data cube via
    the SSA algorithm and return the result.

    Args:
        cube (np.ndarray, ndim=3):
            Data cube which is integrated along the zero-th axis.
    Returns:
        coadded (np.ndarray, ndim=2):
            SSA-integrated frames of the input cube.
    """

    if not isinstance(cube, np.ndarray):
        raise TypeError("specklepy.core.ssa.coadd_frames received cube argument of \
                            type {}, but must be np.ndarray".format(type(cube)))
    if cube.ndim is not 3:
        raise ValueError("specklepy.core.ssa.coadd_frames received cube argument of \
                            dimension {}, but must be 3".format(cube.ndim))

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
    shifts =  shifts.transpose()

    # Shift frames and add to coadded
    coadded = np.zeros(cube[0].shape)
    pad_vectors, ref_pad_vector = alignment.get_pad_vectors(shifts, cube[0].shape, cube[0].shape)
    for index, frame in enumerate(cube):
        coadded += alignment.pad_array(frame, pad_vectors[index], mode='same', reference_image_pad_vector=ref_pad_vector)

    return coadded

