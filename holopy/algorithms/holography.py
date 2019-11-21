import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
import os
import glob
from datetime import datetime
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile
from holopy.core.aperture import Aperture
from holopy.core.apodizer import Apodizer
from holopy.algorithms.psfextraction import PSFExtraction
from holopy.algorithms.sourceextraction import SourceExtraction
from holopy.algorithms.ssa import SSAReconstruction
from holopy.utils.plot import imshow
from holopy.utils.transferfunctions import otf


class HolographicReconstruction(object):

    def __init__(self, params, **kwargs):
        self.params = params
        for key in kwargs:
            self.__setattr__(key, kwargs[key])


    def __call__(self):
        self.execute()


    def execute(self, show=False):
        """Execute the holographic image reconstruction following the algorithm
        outlined in Schoedel et al (2013, Section 3).
        """
        logging.info("Starting holographic reconstruction of {} files...".format(len(self.params.inFiles)))
        self.align_cubes()
        # if not hasattr(self, 'image'):
        self.ssa_reconstruction()
        while True:
            self.find_stars()
            self.select_reference_stars()
            self.extract_psfs()
            self.do_noise_thresholding()
            self.subtract_secondary_sources()
            self.evaluate_object()
            self.apodize_object()
            self.compute_image(autosave=True)

            if show:
                imshow(self.image)

            answer = input("\tDo you want to continue with one more iteration? [yes/no]")
            if answer.lower() in ['n', 'no']:
                break

        return self.image


    def align_cubes(self):
        # self.numberFiles = len(self.params.inFiles)
        # # There is nothing to align, if only one cube is provided.
        # if self.numberFiles == 1:
        #     self.numberFrames = fits.getheader(self.params.inFiles[0])['NAXIS3']
        # else:
        #     self.numberFrames = 0
        #     hdr = fits.getheader(self.params.inFiles[0])
        #     shape = (len(self.params.inFiles), hdr['NAXIS1'], hdr['NAXIS2'])
        #     integrated_cubes = np.zeros(shape)
        #     for index, file in enumerate(self.params.inFiles):
        #         self.numberFrames += fits.getheader(file)['NAXIS3']
        #         integrated_cubes[index] = np.sum(fits.getdata(file), axis=0)
        #     # Implement now how to align the given integrated frames!
        #     raise NotImplementedError("Alignment of multiple cubes is not implemented yet!")
        if len(self.params.inFiles) == 1:
            logging.info("There is nothing to align if only one data cube is provided.")
        else:
            for file in self.params.inFiles:
                integrated = np.sum(fits.getdata(file), axis=0)
                finder = SourceExtraction()
                finder.find_sources(image=integrated, starfinder_fwhm=self.params.starfinderFwhm, noise_threshold=self.params.noiseThreshold,
                    background_subtraction=True, verbose=False)
                if not os.path.isdir(self.params.tmpDir + 'stars/'):
                    os.system('mkdir {}stars/'.format(self.params.tmpDir))
                finder.writeto(self.params.tmpDir + 'stars/' + os.path.basename(file).replace('.fits', '_stars.dat'))



    def ssa_reconstruction(self):
        algorithm = SSAReconstruction()
        self.image = algorithm.execute(self.params.inFiles, outfile=self.params.outFile)
        self.total_flux = np.sum(self.image)


    def find_stars(self, background_subtraction=True):
        """Find sources in the reference image, which may either be a SSA or
        ha preceding holographic reconstruction."""
        finder = SourceExtraction()
        finder.find_sources(image=self.image, starfinder_fwhm=self.params.starfinderFwhm, noise_threshold=self.params.noiseThreshold,
            background_subtraction=background_subtraction, verbose=False)
        finder.writeto(self.params.allStarsFile)


    def select_reference_stars(self):
        """Interactive selection of reference stars"""
        print("\tPlease copy your desired reference stars from the all stars file into the reference star file!")
        input("\tWhen you are done, just hit a key.")


    def extract_psfs(self):
        algorithm = PSFExtraction(self.params)
        algorithm.extract()
        logging.info("Saved the extracted PSFs to the following files:")
        for index, file in enumerate(self.params.psfFiles):
            print('\t', file)


    def do_noise_thresholding(self):
        """Estimate background (np.mean) and noise level (np.std) for every psf
        in a given outer annulus."""
        for file in self.params.psfFiles:
            with fits.open(file, mode='update') as hdulist:
                numberFrames = hdulist[0].header['NAXIS3']
                if not hasattr(self.params, 'psfNoiseMask'):
                    self._generate_noise_mask(hdulist[0].data[0])
                for index in range(numberFrames):
                    reference = np.ma.masked_array(hdulist[0].data[index], mask=self.params.psfNoiseMask)
                    background = np.mean(reference)
                    noise = np.std(reference)
                    update = np.maximum(hdulist[0].data[index] - background - self.params.noiseThreshold * noise, 0.0)
                    update = update / np.sum(update) # Flux sum of order unity
                    hdulist[0].data[index] = update
                    hdulist.flush()


    def _generate_noise_mask(self, frame):
        center = int((frame.shape[0] - 1) / 2)
        radius = center - self.params.noiseReferenceMargin
        tmp = Aperture(center, center, radius, frame, subset_only=False)
        self.params.psfNoiseMask = np.logical_not(tmp.data.mask)


    def subtract_secondary_sources(self):
        pass


    def evaluate_object(self):
        logging.info("Fourier transforming the images...")
        for index, file in enumerate(self.params.inFiles):
            img = fits.getdata(file)
            if index == 0:
                Fimg = fftshift(fft2(img))
            else:
                Fimg = np.concatenate(Fimg, fftshift(fft2(img)))

        logging.info("Fourier transforming the PSFs...")
        for index, file in enumerate(self.params.psfFiles):
            psf = fits.getdata(file)
            if index == 0:
                # Pad the Fpsf cube to have the same xz-extent as Fimg
                print("\tPadding the PSFs")
                print('\tImages', img.shape)
                print('\tPSFs', psf.shape)
                dx = img.shape[1] - psf.shape[1]
                dy = img.shape[2] - psf.shape[2]

                pad_vector = ((0, 0), (int(np.floor(dx/2)), int(np.ceil(dx/2))), (int(np.floor(dy/2)), int(np.ceil(dy/2))))
                print('\tPad_width', pad_vector)
                psf = np.pad(psf, pad_vector)
                try:
                    assert img.shape == psf.shape
                except:
                    raise ValueError("The Fourier transformed images and psfs have different shape, {} and {}. Something went wrong with the padding!".format(Fimg.shape, Fpsf.shape))

                # Initialize Fpsf by the transforming the first cube
                Fpsf = fftshift(fft2(psf))
            else:
                Fpsf = np.concatenate(Fpsf, fftshift(fft2(np.pad(psf, pad_vector))))

        # Clear memory
        del psf
        del img

        # Compute the object
        enumerator = np.mean(np.multiply(Fimg, np.conjugate(Fpsf)), axis=0)
        denominator = np.mean(np.abs(np.square(Fpsf)), axis=0)
        # denominator = np.ma.masked_values(denominator, 0.0)
        self.Fobject  = np.divide(enumerator, denominator)


    def apodize_object(self):
        """Apodize the Fourier object with a apodization function of the users
        choice."""
        # Assert that self.Fobject is square shaped
        try:
            assert self.Fobject.shape[0] == self.Fobject.shape[1]
        except AssertionError:
            raise NotImplementedError("apodization of non-quadratic objects is not implemented yet.")

        logging.info("Apodizing the object...")
        self.apodizer = Apodizer(self.params.apodizationType, size=self.Fobject.shape[0], radius=self.params.apodizationWidth)
        self.Fobject = self.apodizer.apodize(self.Fobject)


    def compute_image(self, autosave=True):
        """Compute the final image from the apodized object."""
        self.image = ifft2(self.Fobject)
        self.image = np.abs(self.image)
        image_scale = self.total_flux / np.sum(self.image)
        self.image = np.multiply(self.image, image_scale) # assert flux conservation
        if autosave:
            self.params.outFile.data = self.image
        return self.image
