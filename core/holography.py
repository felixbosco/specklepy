import numpy as np
from scipy import signal
import os
import glob
from datetime import datetime
from astropy.io import fits
import logging
from logging.config import fileConfig
fileConfig('./config/logging.cfg')

from lib.reconstructor import Reconstructor
from lib.visual import imshow


class HolographicReconstructor(Reconstructor):

    def __init__(self, **kwargs):
        super(HolographicReconstructor, self).__init__(**kwargs)


    def reconstruct(self, **kwargs):
        if not 'object_file' in kwargs:
            self.estimate_psfs()
            object = self.estimate_object()
        else:
            try:
                real, imag = fits.getdata(kwargs['object_file'])
            except:
                real, imag = fits.getdata(self.tmp + kwargs['object_file'])
            object = np.empty(real.shape, dtype=complex)
            object.real = real
            object.imag = imag

        apodized = self.apodization(object)
        reconstruction = np.fft.fftshift(np.abs(np.fft.ifft2(apodized)))

        # save the result to a file in self.output
        hdu = fits.PrimaryHDU(reconstruction)
        hdu.header['OBJECT'] = 'Holopy holographic reconstruction'
        for index, file in enumerate(glob.glob(self.input + self.cube_file)):
            hdu.header['HIERARCH HOLOPY FILE {}'.format(index)] = os.path.basename(file)
        hdu.header['DATE'] = str(datetime.now())
        hdulist = fits.HDUList([hdu])
        logging.info("Saving result to: {}".format(self.output + 'holo.fits'))
        hdulist.writeto(self.output + 'holo.fits', overwrite=True)
        logging.info("Saving successful!")

        return reconstruction


    def estimate_psfs(self):
        pass


    def estimate_object(self):
        # creating file lists
        image_file_list = glob.glob(self.input + self.cube_file)
        psf_file_list = [self.tmp + os.path.basename(file).replace('.fits', '_psf.fits') for file in  image_file_list]

        # check input files
        with fits.open(image_file_list[0]) as hdulist:
            image_shape = hdulist[0].data.shape
        with fits.open(psf_file_list[0]) as hdulist:
            psf_shape = hdulist[0].data.shape
        logging.info("Computing the object reconstruction with image and PSF data of shapes:")
        logging.info(image_shape)
        logging.info(psf_shape)
        if image_shape[0] != psf_shape[0]:
            raise ValueError("Image and PSF input have different number of frames!")

        # estimating required zero-padding
        logging.info("Padding the PSF data to have the same shape as the image data...")
        pad_width = self._get_pad_width(image_shape, psf_shape)

        # computing enumerator and denominator separatly.
        # both have the same frame number, therefore the average is equivalent to
        # a simple sum...
        enumerator = np.zeros((image_shape[1], image_shape[2]), dtype=complex)
        denominator = np.zeros((image_shape[1], image_shape[2]), dtype=complex)
        logging.info("Computing the enumerator...")
        for file_index in range(len(image_file_list)):
            image_input = fits.getdata(image_file_list[file_index])
            psf_input = fits.getdata(psf_file_list[file_index])
            for index, frame in enumerate(image_input):
                print("File {:2} Frame {:4} of {:4}...".format(file_index + 1, index + 1, image_input.shape[0]), end='\r')
                image = np.fft.fft2(frame)
                psf = np.fft.fft2(np.pad(psf_input[index], pad_width=pad_width, mode='constant'))
                enumerator += np.multiply(image, psf)
            logging.info("Computing the denominator...")
            for index, frame in enumerate(psf_input):
                print("File {:2} Frame {:4} of {:4}...".format(file_index + 1, index + 1, psf_input.shape[0]), end='\r')
                psf = np.fft.fft2(np.pad(frame, pad_width=pad_width, mode='constant'))
                denominator += np.abs(np.square(psf))

        if False:
            enumerator = np.fft.fftshift(enumerator)
            denominator = np.fft.fftshift(denominator)

        if False:
            imshow([enumerator.real, enumerator.imag, denominator.real, denominator.imag],
                        title=['enumerator.real', 'enumerator.imag', 'denominator.real', 'denominator.imag'],
                        figsize=(15, 5))

        object = np.divide(enumerator, denominator)

        if hasattr(self, 'save_tmp') and self.save_tmp:
            fits.writeto(self.tmp + 'object.fits', data=np.array([object.real, object.imag]), overwrite=True)

        return object


    def apodization(self, object):
        """
        This function computes the apodization with a telescope optical transfer
        function (OTF), obtained from the Fourier transform of a given telescope
        point spread function(PSF) - Step (x) of the holography algorithm.
        """
        logging.info("Reading PSF from file {}...".format(self.psf_file))
        out_psf = fits.getdata(self.input + self.psf_file)
        logging.info("Apodizing the object with the given PSF in the Fourier domain...")
        otf = np.fft.fftshift(np.fft.fft2(out_psf))
        return signal.convolve2d(object, otf, mode='same')
        return object


    def _get_pad_width(self, shape_large, shape_small):
        shape_diff = [shape_large[dim] - shape_small[dim] for dim in range(1, 3)]
        pad_width = [(int(shape_diff[dim] / 2), int(np.ceil(shape_diff[dim] / 2))) for dim in range(2)]
        return pad_width
