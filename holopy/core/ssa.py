import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile
from holopy.core.alignment import compute_shifts


def ssa(file_list, mode='same', reference_file=None, reference_file_index=0, outfile=None, **kwargs):
    """Compute the SSA reconstruction of a list of files.

    Long description...

    Args:
        file_list (list):
        mode (str):
        reference_file (str, optional): Path to a reference file, relative to
            which the shifts are computed. If not provided, the reference file
            index is used. See holopy.core.aligment.compute_shifts for details.
        reference_file_index (str, optional): Index of the file in the file
            list, relative to which the shifts are cpmputed. See
            holopy.core.aligment.compute_shifts for details. Default is 0.
        outfile (holopy.io.outfile, optional): Object to write the result to,
            if provided.
    """
    logging.info("Starting SSA reconstruction...")

    file_shifts, image_shape = compute_shifts(file_list, reference_file=reference_file, reference_file_index=reference_file_index, return_image_shape=True, lazy_mode=True)
    reconstruction = np.zeros(image_shape)

    for index, file in enumerate(file_list):
        cube = fits.getdata(file)
        reconstruction = reconstruction + shift_array(coadd_frames(cube), shift=file_shifts[index])

    logging.info("Reconstruction finished...")

    # Save the result to an Outfile
    if outfile is not None:
        outfile.data = reconstruction

    return reconstruction


def coadd_frames(cube, mode='same'):
    """
    Compute the simple shift-and-add (SSA) reconstruction via the SSA algorithm
    of a fits cube and return the result.
    """

    # Compute shifts
    peak_indizes = np.zeros((cube.shape[0], 2), dtype=int)
    for index, frame in enumerate(cube):
        peak_indizes[index] = np.array(np.unravel_index(np.argmax(frame, axis=None), frame.shape), dtype=int)
    shifts = compute_shifts_from_indizes(peak_indizes)


    # Shift frames and add to out
    out = np.zeros(cube[0].shape)
    for index, frame in enumerate(cube):
        out += shift_array(frame, shifts[index])

    return out



def compute_shifts_from_indizes(indizes):
    indizes = indizes.transpose()
    xmean, ymean = np.mean(np.array(indizes), axis=1)
    xmean = int(xmean)
    ymean = int(ymean)
    shifts = np.array([indizes[0] - xmean, indizes[1] - ymean])
    return shifts.transpose()



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
