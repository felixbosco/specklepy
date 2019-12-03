import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from specklepy.logging import logging
from specklepy.io.outfile import Outfile
from specklepy.core import alignment


def ssa(files, mode='same', reference_file=0, outfile=None, tmp_dir=None, lazy_mode=True, debug=False, **kwargs):
    """Compute the SSA reconstruction of a list of files.

    Long description...

    Args:
        files (list):
        mode (str):
        reference_file (str, int, optional):
            Path to a reference file or index of the file in files, relative to
            which the shifts are computed. See
            specklepy.core.aligment.get_shifts for details. Default is 0.
        outfile (specklepy.io.outfile, optional): Object to write the result to,
            if provided.
        debug (bool, optional): Set to True to inspect intermediate results.
            Default is False.

    Returns:
        reconstruction (np.ndarray): The image reconstruction.
    """

    logging.info("Starting SSA reconstruction...")
    if not isinstance(files, list):
        files = [files]

    if isinstance(reference_file, int):
        reference_file = files[reference_file]
    elif not isinstance(reference_file, str):
        raise TypeError("The function get_shifts received reference_file argument of type {}, but needs be int or str, i.e. a file name.".format(type(reference_file)))

    if outfile is not None and not isinstance(outfile, Outfile):
        if isinstance(outfile, str):
            outfile = Outfile(files=files, filename=outfile, cards={"RECONSTRUCTION": "SSA"})
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

    Returns:
        coadded (np.ndarray, ndim=2): SSA-integrated frames of the input cube.
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
    # shifts = alignment_from_indizes(peak_indizes)

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



def create_pad_vector_entry(shift_entry):
    if shift_entry <= 0 :
        return (np.abs(shift_entry), 0)
    else:
        return (0, shift_entry)



def create_pad_vector(shift):
    return (create_pad_vector_entry(shift[0]), create_pad_vector_entry(shift[1]))



def shift_array(array, shift):
    shape = array.shape
    pad_array = array[max(0, shift[0]) : shape[0]+min(0, shift[0]) , max(0, shift[1]) : shape[1]+min(0, shift[1])]
    pad_vector = create_pad_vector(shift)
    return np.pad(pad_array, pad_vector, mode='constant')
