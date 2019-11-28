import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile
from holopy.core.alignment import compute_shifts


class SSAReconstruction(object):

    def __init__(self, **kwargs):
        # super(SSAReconstruction, self).__init__(**kwargs)
        for key in kwargs:
            self.__setattr__(key, kwargs[key])


    def __call__(self, *args, **kwargs):
        self.execute(*args, **kwargs)


    def execute(self, file_list, outfile=None, **kwargs):
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
        logging.info("Starting {}...".format(self.__class__.__name__))

        for index, file in enumerate(file_list):
            cube = fits.getdata(file)

            if index == 0:
                reconstruction = self._ssa(cube)
            else:
                reconstruction = self._align_reconstructions(reconstruction, self._ssa(cube), shift=None)#file_shifts[index])

        logging.info("Reconstruction finished...")

        # Save the result to an Outfile
        if outfile is not None:
            outfile.data = reconstruction

        return reconstruction


    def _ssa(self, cube):
        """
        Compute the simple shift-and-add (SSA) reconstruction via the SSA algorithm
        of a fits cube and return or write the result to a file.
        """

        # Initialize
        indizes = []
        out = np.zeros(cube[0].shape)

        # Compute shifts
        for index, frame in enumerate(cube):
            print('Estimating the shift in frame {:4}...'.format(index + 1), end='\r')
            indizes.append(np.array(np.unravel_index(np.argmax(frame, axis=None), frame.shape)))
        shifts = self._compute_shifts(indizes)
        print('Estimated all shifts                           ', end='\r')

        # Shift frames and add to out
        for index, frame in enumerate(cube):
            print('Adding the frame {:4}...'.format(index + 1), end='\r')
            out += self._shift_array(frame, shifts[index])

        return out


    def _compute_shifts(self, indizes):
        xmean, ymean = np.mean(np.array(indizes), axis=0)
        xmean = int(xmean)
        ymean = int(ymean)
        shifts = []
        for i in indizes:
            shifts.append((i[0]-xmean, i[1]-ymean))
        return shifts


    def _create_pad_vector_entry(self, shift_entry):
        if shift_entry <= 0 :
            return (np.abs(shift_entry), 0)
        else:
            return (0, shift_entry)


    def _create_pad_vector(self, shift):
        return (self._create_pad_vector_entry(shift[0]), self._create_pad_vector_entry(shift[1]))


    def _shift_array(self, array, shift):
        shape = array.shape
        pad_array = array[max(0, shift[0]) : shape[0]+min(0, shift[0]) , max(0, shift[1]) : shape[1]+min(0, shift[1])]
        pad_vector = self._create_pad_vector(shift)
        return np.pad(pad_array, pad_vector, mode='constant')


    def _align_reconstructions(self, previous, current, shift=None):
        if shift is None:
            previous_x0, previous_y0 = np.unravel_index(np.argmax(previous, axis=None), previous.shape)
            current_x0, current_y0 = np.unravel_index(np.argmax(current, axis=None), current.shape)
            shift = (current_x0 - previous_x0, current_y0 - previous_y0)
        return previous + self._shift_array(current, shift)
