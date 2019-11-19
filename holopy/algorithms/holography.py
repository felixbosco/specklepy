import numpy as np
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile
from holopy.core.apodizer import Apodizer
from holopy.algorithms.psfextraction import PSFExtraction


class HolographicReconstruction(object):

    def __init__(self, params, **kwargs):
        self.params = params
        for key in kwargs:
            self.__setattr__(key, kwargs[key])

    def __call__(self):
        self.execute()


    def execute(self):
        """Execute the holographic image reconstruction following the algorithm
        outlined in Schoedel et al (2013, Section 3).
        """
        logging.info("Starting holographic reconstruction of {} files...".format(len(self.params.inFiles)))
        self.align_cubes()
        self.find_stars()
        self.select_reference_stars()
        self.extract_psfs()
        self.do_noise_thresholding()
        self.subtract_secondary_sources()
        self.evaluate_object()
        self.apodize_object()
        self.compute_image()


    def align_cubes(self):
        pass


    def find_stars(self):
        """Find sources in the reference image, which may either be a SSA or
        ha preceding holographic reconstruction."""
        pass


    def select_reference_stars(self):
        pass


    def extract_psfs(self):
        algorithm = PSFExtraction(self.params)
        algorithm.extract()
        logging.info("Saved the extracted PSFs to the following files:")
        for file in self.params.psfFiles:
            print(file)


    def do_noise_thresholding(self):
        pass


    def subtract_secondary_sources(self):
        pass


    def evaluate_object(self, images, psfs, apodizer):
        try:
            assert images.shape == psfs.shape
        except:
            raise NotImplementedError('Padding of the PSF estimate for obtaining the proper shape is not implemented yet')

        Fimg = np.fft.fft2(images)
        Fpsf = np.fft.fft2(psfs)

        Fobj = np.divide(np.mean(np.multiply(Fimg, np.conjugate(Fpsf)), axis=0), np.mean(np.abs(np.square(Fpsf)), axis=0))
        self.Fobject = Fobj


    def apodize_object(self):
        """Apodize the Fourier object with a apodization function of the users
        choice."""
        # Assert that self.Fobject is square shaped
        try:
            assert self.Fobject.shape[0] == self.Fobject.shape[1]
        except AssertionError:
            raise NotImplementedError("apodization of non-quadratic objects is not implemented yet.")

        self.apodizer = Apodizer(params.apodizationType, size=self.Fobject.shape[0], radius=params.apodizationWidth)
        self.Fobject = self.apodizer.apodize(self.Fobject)


    def compute_image(self, autosave=True):
        """Compute the final image from the apodized object."""
        self.image = np.fft.ifft2(self.Fobject)
        if autosave:
            with fits.open(self.params.outFile, mode='update') as hdulist:
                # IDEA implement that preceeding reconstructions are saved into the extensions
                hdulist[0].data = self.images
                hdulist.flush()
        return self.image
