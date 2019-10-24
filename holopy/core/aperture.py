import numpy as np
from copy import copy
from astropy.io import fits
from datetime import datetime

from holopy.logging import logging
from holopy.utils import transferfunctions as tf


class Aperture(object):

    def __init__(self, x0, y0, radius, data, mask='circular', subset_only=True):
        self.x0 = x0
        self.y0 = y0
        self.radius = radius

        # Handling data input
        if isinstance(data, str):
            logging.info("Aperture argument data '{}' is interpreted as file name.".format(data))
            data = fits.getdata(data)
        if not (data.ndim == 2 or data.ndim == 3):
            raise ValueError("Data input of Aperture class must be of dimension 2 or 3, but was provided as data.ndim={}.".format(data.ndim))
        self.data = copy(data)

        # Interprete mask argument
        self.mask = mask
        if self.mask is None:
            pass
        elif self.mask is 'circular':
            self.data = np.ma.masked_array(self.data, mask=self.make_mask())
        else:
            raise ValueError("Mask type '{}' of Aperture instance not understood.".format(mask))

        # Remove the masked margins if requested
        self.subset_only = subset_only
        if subset_only:
            if data.ndim == 2:
                self.data = copy(self.data[self.x0 - self.radius : self.x0 + self.radius + 1, self.y0 - self.radius : self.y0 + self.radius + 1])
            elif data.ndim == 3:
                self.data = copy(self.data[:, self.x0 - self.radius : self.x0 + self.radius + 1, self.y0 - self.radius : self.y0 + self.radius + 1])


    @property
    def width(self):
        return self.radius * 2 + 1


    @property
    def index(self):
        return (self.x0, self.y0)


    def __call__(self):
        return self.data


    def make_distance_map(self):
        xx, yy = np.mgrid[:self.data.shape[-2], :self.data.shape[-1]]
        return np.sqrt(np.square(xx - self.x0) + np.square(yy - self.y0))


    def make_mask(self):
        distance_map = self.make_distance_map()
        mask2D = np.ma.masked_greater(distance_map, self.radius).mask
        if self.data.ndim == 2:
            return mask2D
        elif self.data.ndim == 3:
            mask_3D = np.expand_dims(mask2D, axis=0)
            return np.repeat(mask3D, repeats=self.data.shape[0], axis=0)


    def get_aperture_peak(self):
        if self.data.ndim == 3:
            logging.info("Aperture data are 3D and therefore integrated before the aperture is recentered.")
            # Integrate the data cube
            tmp = np.sum(self.data, axis=0)
            return np.unravel_index(np.argmax(tmp, axis=None), tmp.shape)
        else:
            return np.unravel_index(np.argmax(self.data, axis=None), self.data.shape)


    def recenter_aperture(self):
        if self.subset_only:
            logging.warning("Recentering the aperture on the peak within the guess aperture is working only if subset_only is set to False.")
        else:
            peak = self.get_aperture_peak()
            self = self.__init__(*peak, self.radius, data=self.data.getdata(), mask=self.mask, subset_only=self.subset_only)


    def initialize_Fourier_file(self, infile, Fourier_file):
        self.infile = infile
        self.Fourier_file = Fourier_file
        logging.info("Initializing Fourier file {}".format(self.Fourier_file))
        header = fits.getheader(self.infile)
        header.set('HIERARCH HOLOPY TYPE', 'Fourier transform of an aperture')
        header.set('HIERARCH HOLOPY ORIGIN', self.infile)
        header.set('HIERARCH HOLOPY APERTURE INDEX', str(self.index))
        header.set('HIERARCH HOLOPY APERTURE RADIUS', self.radius)
        header.set('UPDATED', str(datetime.now()))
        data = np.zeros(self.data.shape)
        fits.writeto(self.Fourier_file, data=data, header=header, overwrite=True)
        logging.info("Initialized {}".format(self.Fourier_file))


    def powerspec_to_file(self, infile=None, Fourier_file=None):
        if not hasattr(self, 'Fourier_file'):
            self.initialize_Fourier_file(infile, Fourier_file)

        with fits.open(self.Fourier_file, mode='update') as hdulist:
            for index, frame in enumerate(self.data):
                print("\rFourier transforming frame {}/{}".format(index+1, self.data.shape[0]), end='')
                hdulist[0].data[index] = tf.powerspec(frame)
                hdulist.flush()
            print()
        logging.info("Computed the Fourier transform of every frame and saved them to {}".format(self.Fourier_file))


    def powerspec(self):
        self.Fourier_data = np.zeros(self.data.shape)
        for index, frame in enumerate(self.data):
            print("\rFourier transforming frame {}/{}".format(index+1, self.data.shape[0]), end='')
            self.Fourier_data[index] = tf.powerspec(frame)
        print()
        logging.info("Computed the Fourier transform of every frame.")


    def get_encircled_energy(self, saveto=None):
        if self.data.ndim == 3:
            logging.info("Aperture data are 3D and integrated before the encircled energy is estimated.")
            # Integrate the data cube
            self.tmp = np.sum(self.data, axis=0)
        else:
            self.tmp = copy(self.data)

        # logging.warning('encircled_energy() is not implemented yet.')

        rdata = np.arange(0, self.radius, 1)
        out = np.zeros((self.radius))
        distance_map = self.make_distance_map()
        for index, subset_radius in enumerate(rdata):
            out[index] = np.sum(self.tmp[np.where(distance_map <= subset_radius)])

        if saveto is not None:
            # Save power spectra to outfile
            caption = "radius_(pix) encircled_energy(data_unit)"
            data = np.concatenate(([rdata], [out]), axis=0).transpose()
            np.savetxt(saveto, data, header=caption)

        del self.tmp
        return rdata, out
