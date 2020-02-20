import numpy as np
import warnings
from copy import copy
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logger
from specklepy.utils import transferfunctions as tf


class Aperture(object):

    def __init__(self, *args, data, vars=None, mask='circular', crop=True, verbose=True):
        """Instantiate Aperture object.

        Long description...

        Args:
            *args (int, float:
                Can be either two or three arguments, where the last is always
                interpreted as the aperture radius, which must be int type. The
                other arguments can be either a coordinate tuple or two
                individual coordinates.
            data (np.ndarray, str):
                2D or 3D np.ndarray that the aperture shall be extracted of. If
                provided as str type, this is assumed to be a path and the
                objects tries to read the fits file.
            vars (np.ndarray, optional):
                Variance of the data. If not provided, the variance will be
                estimated along the time axis of a cube. In the future this
                might read in from a file extension. Default is None.
            mask (str, optional):
                Mode that is describing, how the aperture is masked. Can be
                'circular' or 'rectangular'. If 'circular', then it creates a
                circular mask and the data become a np.ma.masked_array. Default
                is 'circular'.
            crop (bool, optional):
                If set to True, then the object only stores a copy of the data
                with radius around the center. Otherwise all the data beyond
                the limits of the aperture are masked. Default is True.
            verbose (bool, optional):
                Set to True for retrieving more information. Default is True.
        """

        # Interprete the args
        if len(args) == 2 and (isinstance(args[0], tuple) or isinstance(args[0], list)):
            x0 = args[0][0]
            y0 = args[0][1]
            radius = args[1]
        elif len(args) == 3:
            x0 = args[0]
            y0 = args[1]
            radius = args[2]
        else:
            raise ValueError("Aperture expects either 2 or 3 args, but {} have been provided!".format(len(args)))

        # Integer checks and handling non-integer input
        if isinstance(x0, int) and isinstance(y0, int):
            self.x0 = x0
            self.y0 = y0
            self.xoffset = 0
            self.yoffset = 0
        else:
            x0 = float(x0)
            y0 = float(y0)
            self.x0 = np.rint(x0).astype(int)
            self.y0 = np.rint(y0).astype(int)
            self.xoffset = x0 - self.x0
            self.yoffset = y0 - self.y0


        if isinstance(radius, int):
            self.radius = radius
        else:
            raise ValueError('Aperture radius must be given as integer! (Was given as {})'.format(radius))

        # Handling data input
        if isinstance(data, str):
            if verbose:
                logger.info("Aperture argument data '{}' is interpreted as file name.".format(data))
            data = fits.getdata(data)
        if not (data.ndim == 2 or data.ndim == 3):
            raise ValueError("Data input of Aperture class must be of dimension 2 or 3, but was provided as data.ndim={}.".format(data.ndim))
        self.data = copy(data)

        # Remove the margins if requested
        self.cropped = False
        if crop:
            self.crop()

        # Create a mask
        mask = self.make_mask(mode=mask)
        self.data = np.ma.masked_array(self.data, mask=mask)

        if vars is None:
            if self.data.ndim == 3:
                self.vars = np.var(self.data, axis=0)
        else:
            self.vars = None



    @property
    def width(self):
        return self.radius * 2 + 1

    @property
    def index(self):
        return (self.x0, self.y0)

    @property
    def offset(self):
        return (self.xoffset, self.yoffset)

    @property
    def shape(self):
        return self.data.shape

    @property
    def nframes(self):
        return self.data.shape[0]

    def __getitem__(self, index):
        return self.data[index]



    def make_distance_map(self):
        if self.cropped:
            center = (self.radius, self.radius)
        else:
            center = self.index
        xx, yy = np.mgrid[:self.data.shape[-2], :self.data.shape[-1]]
        return np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))



    def make_mask(self, mode='circular'):
        """Create a circular or rectangular mask."""

        if not isinstance(mode, str):
            raise TypeError("Aperture received mask argument of type {}, but needs to be str type!".format(type(mask)))

        if mode == 'circular':
            distance_map = self.make_distance_map()
            mask2D = np.ma.masked_greater(distance_map, self.radius).mask
            if self.data.ndim == 2:
                return mask2D
            elif self.data.ndim == 3:
                mask3D = np.expand_dims(mask2D, axis=0)
                return np.repeat(mask3D, repeats=self.data.shape[0], axis=0)
        elif mode == 'rectangular':
            if self.cropped:
                return np.zeros(self.data.shape, dtype=bool)
            else:
                mask = np.ones(self.data.shape, dtype=bool)
                if mask.ndim == 2:
                    mask[self.x0 - self.radius : self.y0 + self.radius + 1, self.y0 - self.radius : self.y0 + self.radius +1] = 0
                elif mask.ndim == 3:
                    mask[ : , self.x0 - self.radius : self.y0 + self.radius + 1, self.y0 - self.radius : self.y0 + self.radius +1] = 0
                return mask
        else:
            raise ValueError("Aperture received mask argument {}, but needs to be either 'circular' or 'rectangular'!".format(mask))



    def crop(self):
        if self.cropped:
            logger.info("Margins are removed already from aperture instance.")
        else:
            if self.data.ndim == 2:
                self.data = copy(self.data[self.x0 - self.radius : self.x0 + self.radius + 1, self.y0 - self.radius : self.y0 + self.radius + 1])
            elif self.data.ndim == 3:
                self.data = copy(self.data[:, self.x0 - self.radius : self.x0 + self.radius + 1, self.y0 - self.radius : self.y0 + self.radius + 1])
            self.cropped = True

    def remove_margins(self):
        self.crop()



    def get_integrated(self):
        """Returns a 2-dimensional represntation of the aperture and integrates
        along the time axis if necessary.
        """
        if self.data.ndim == 3:
            return np.sum(self.data, axis=0)
        else:
            return copy(self.data)



    def get_aperture_peak(self):
        """Returns the coordinates of the emission peak in the aperture."""
        tmp = self.get_integrated()
        return np.unravel_index(np.argmax(tmp, axis=None), tmp.shape)



    def initialize_Fourier_file(self, infile, Fourier_file):
        self.infile = infile
        self.Fourier_file = Fourier_file
        logger.info("Initializing Fourier file {}".format(self.Fourier_file))
        header = fits.getheader(self.infile)
        header.set('HIERARCH specklepy TYPE', 'Fourier transform of an aperture')
        header.set('HIERARCH specklepy ORIGIN', self.infile)
        header.set('HIERARCH specklepy APERTURE INDEX', str(self.index))
        header.set('HIERARCH specklepy APERTURE RADIUS', self.radius)
        header.set('UPDATED', str(datetime.now()))
        data = np.zeros(self.data.shape)
        fits.writeto(self.Fourier_file, data=data, header=header, overwrite=True)
        logger.info("Initialized {}".format(self.Fourier_file))



    def powerspec_to_file(self, infile=None, Fourier_file=None):
        if not hasattr(self, 'Fourier_file'):
            self.initialize_Fourier_file(infile, Fourier_file)

        with fits.open(self.Fourier_file, mode='update') as hdulist:
            for index, frame in enumerate(self.data):
                print("\rFourier transforming frame {}/{}".format(index+1, self.data.shape[0]), end='')
                hdulist[0].data[index] = tf.powerspec(frame)
                hdulist.flush()
            print()
        logger.info("Computed the Fourier transform of every frame and saved them to {}".format(self.Fourier_file))



    def powerspec(self):
        self.Fourier_data = np.zeros(self.data.shape)
        for index, frame in enumerate(self.data):
            print("\rFourier transforming frame {}/{}".format(index+1, self.data.shape[0]), end='')
            self.Fourier_data[index] = tf.powerspec(frame)
        print()
        logger.info("Computed the Fourier transform of every frame.")



    def get_encircled_energy(self, saveto=None):
        """Extracts the encircled energy from an aperture as a function of
        radius."""

        tmp = self.get_integrated()

        # Initialize variables
        rdata = np.arange(0, self.radius, 1)
        ydata = np.zeros(rdata.shape)
        distance_map = self.make_distance_map()

        # Iterate over aperture radii
        for index, subset_radius in enumerate(rdata):
            mask = np.ma.masked_greater(distance_map, subset_radius).mask
            ydata[index] = np.sum(np.ma.masked_array(tmp, mask=mask))

        # Save results to file
        if saveto is not None:
            header = "radius_(pix) encircled_energy(data_unit)"
            data = np.concatenate(([rdata], [ydata]), axis=0).transpose()
            np.savetxt(saveto, data, header=header)
            logger.info("Saved encircled energy data to file {}".format(saveto))

        return rdata, ydata



    def get_psf_profile(self):
        """Computes the radial average of the PSF in the aperture."""

        tmp = self.get_integrated()

        radius_map = self.make_distance_map()
        rdata = np.unique(radius_map)
        ydata = np.zeros(rdata.shape)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for index, radius in enumerate(rdata):
                ydata[index] = np.mean(tmp[np.where(radius_map == radius)])

        return rdata, ydata
