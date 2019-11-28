import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile
from holopy.core.alignment import compute_shifts


def ssa(files, mode='same', reference_file=None, reference_file_index=0, outfile=None, tmp_dir=None, lazy_mode=True, debug=False, **kwargs):
    """Compute the SSA reconstruction of a list of files.

    Long description...

    Args:
        files (list):
        mode (str):
        reference_file (str, optional): Path to a reference file, relative to
            which the shifts are computed. If not provided, the reference file
            index is used. See holopy.core.aligment.compute_shifts for details.
        reference_file_index (str, optional): Index of the file in the file
            list, relative to which the shifts are cpmputed. See
            holopy.core.aligment.compute_shifts for details. Default is 0.
        outfile (holopy.io.outfile, optional): Object to write the result to,
            if provided.
        debug (bool, optional): Set to True to inspect intermediate results.
            Default is False.

    Returns:
        reconstruction (np.ndarray): The image reconstruction.
    """

    logging.info("Starting SSA reconstruction...")
    if not isinstance(files, list):
        files = [files]
    if outfile is not None and not isinstance(outfile, Outfile):
        if isinstance(outfile, str):
            outfile = Outfile(files=files, filename=outfile, cards={"RECONSTRUCTION": "SSA"})
        else:
            raise TypeError("holopy.core.ssa.ssa received outfile argument of wrong type <{}>!".format(type(outfile)))


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
        file_shifts, image_shape = compute_shifts(tmp_files, reference_file=reference_file, reference_file_index=reference_file_index, return_image_shape=True, lazy_mode=True)
        reconstruction = np.zeros(image_shape)
        for index, file in enumerate(tmp_files):
            tmp_image = fits.getdata(file)
            reconstruction = reconstruction + shift_array(tmp_image, shift=file_shifts[index])

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
    # shifts = compute_shifts_from_indizes(peak_indizes)

    # Compute shifts from indizes
    peak_indizes = peak_indizes.transpose()
    xmean, ymean = np.mean(np.array(peak_indizes), axis=1)
    xmean = int(xmean)
    ymean = int(ymean)
    shifts = np.array([peak_indizes[0] - xmean, peak_indizes[1] - ymean])
    shifts =  shifts.transpose()

    # Shift frames and add to out
    out = np.zeros(cube[0].shape)
    for index, frame in enumerate(cube):
        out += shift_array(frame, shifts[index])

    return out



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
