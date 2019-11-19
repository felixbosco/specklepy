import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile


class HolographicReconstruction(object):

    def __init__(self, params, **kwargs):
        self.params = params
        for key in kwargs:
            self.__setattr__(key, kwargs[key])


    def __call__(self, file_list, outfile=None):
        logging.info("Starting {}...".format(self.__class__.__name__))
        # file_list = glob.glob(self.input + self.cube_file)

        for index, file in enumerate(file_list):
            cube = fits.getdata(file)
            if index == 0:
                reconstruction = self._ssa(cube)
            else:
                reconstruction = self._align_reconstructions(reconstruction, self._ssa(cube))

        logging.info("Reconstruction finished...")

        # Save the result to an Outfile
        if outfile is not None:
            outfile.data = reconstruction

        return reconstruction


    def align_cubes(self):
        pass


    def find_stars(self):
        """Find sources in the reference image, which may either be a SSA or
        ha preceding holographic reconstruction."""
        pass


    def select_reference_stars(self):
        pass


    def estimate_psfs(self):
        pass


    def noise_thresholding(self):
        pass


    def secondary_source_subtraction(self):
        pass


    def evaluate_object(self, images, psfs, apodizer):
        try:
            assert images.shape == psfs.shape
        except:
            raise NotImplementedError('Padding of the PSF estimate for obtaining the proper shape is not implemented yet')

        Fimg = np.fft.fft2(images)
        Fpsf = np.fft.fft2(psfs)

        Fobj = np.divide(np.mean(np.multiply(Fimg, np.conjugate(psfs), np.mean(..., axis=0)), axis=0))
        self.Fobject = Fobj


    def apodize(self):
        """Apodize the Fourier object with a apodization function of the users
        choice."""
        pass


    def compute_image(self):
        """Compute the final image from the apodized object."""
        self.image = np.fft.fft2(self.Fobject)
        return self.image
