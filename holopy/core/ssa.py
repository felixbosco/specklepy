import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile
from holopy.core.reconstructor import Reconstructor


class SSAReconstructor(Reconstructor):

    def __init__(self, **kwargs):
        # super(SSAReconstructor, self).__init__(**kwargs)
        for key in kwargs:
            self.__setattr__(key, kwargs[key])


    def __call__(self, file_list):
        logging.info("Starting reconstruction with {}...".format(self.__class__.__name__))
        # file_list = glob.glob(self.input + self.cube_file)

        for index, file in enumerate(file_list):
            cube = fits.getdata(file)
            if index == 0:
                reconstruction = np.zeros(cube[0].shape)

            current_reconstruction = self._ssa(cube)
            reconstruction += current_reconstruction

            save_to = os.path.basename(file)
            save_to = save_to.replace('.fits', '_ssa.fits')
            if hasattr(self, 'save_tmp') and self.save_tmp:
                fits.writeto(self.tmp + save_to, data=current_reconstruction)
        logging.info("Reconstruction finished...")

        # save the result to a file in self.output
        self.outfile.data = reconstruction
        # hdu = fits.PrimaryHDU(reconstruction)
        # hdu.header['OBJECT'] = 'Holopy SSA reconstruction'
        # for index, file in enumerate(file_list):
        #     hdu.header['HIERARCH HOLOPY FILE {}'.format(index)] = os.path.basename(file)
        # hdu.header['DATE'] = str(datetime.now())
        # hdulist = fits.HDUList([hdu])
        # logging.info("Saving result to: {}".format(self.output + 'ssa.fits'))
        # hdulist.writeto(self.output + 'ssa.fits', overwrite=True)
        # logging.info("Saving successful!")

        return reconstruction


    def _ssa(self, cube):
        """
        Compute the simple shift-and-add (SSA) reconstruction via the SSA algorithm
        of a fits cube and return or write the result to a file.
        """

        # initialize
        indizes = []
        out = np.zeros(cube[0].shape)

        # compute shifts
        for index, frame in enumerate(cube):
            print('Estimating the shift in frame {:4}.'.format(index + 1), end='\r')
            indizes.append(np.array(np.unravel_index(np.argmax(frame, axis=None), frame.shape)))
        shifts = self._compute_shifts(indizes)

        # shift frames and add to out
        for index, frame in enumerate(cube):
            print('Adding the frame {:4}.'.format(index + 1), end='\r')
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
