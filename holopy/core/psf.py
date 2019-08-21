import os
import sys
import numpy as np
import glob
from astropy.io import fits

from holopy.logging import logging

# from lib.rectangle import Rectangle
# from lib.refsource import RefSource


class PsfEstimator(object):

    def __init__(self, **kwargs):
        # read parameter file
        for key in ['parameterFileName', 'parameter_file_name']:
            if key in kwargs:
                self.read_parameter_file(kwargs[key])
        # save kwargs to __dict__
        for key in kwargs:
            setattr(self, key, kwargs[key])
        # apply defaults to unset but required parameters
        self.apply_defaults()
        self.create_directories()
        self.initialize_psf_files()


    def __str__(self):
        str = "PsfEstimator"
        str += '\n' + len(str) * '-' + '\n'
        for key in self.__dict__:
            str += "\n{} = {}".format(key, self.__dict__[key])
        return str


    def __call__(self, **kwargs):
        logging.info('Estimating PSFs in the input files...')
        self.estimate_psf(**kwargs)


    def save(self):
        with fits.open(self.out_file_name, mode='update') as hdulist:
            hdulist[0].data[self.current_frame] = self.current_psf
        hdulist.flush()


    def read_parameter_file(self, file_name, delimiter='=', comment='#'):
        logging.info('Reading in parameter file...')
        self.parameter_file_name = file_name
        with open(file_name) as f:
            content = f.readlines()
            for line in content:
                # remove breaks
                line = line.replace('\n', '')
                # remove spaces
                line = line.replace(' ', '')
                # remove comments
                line = line.split(comment)[0]
                if len(line) == 0 :
                    continue
                # assign key and value to self.__dict__
                key, value = line.split(delimiter)
                setattr(self, key, value)


    def read_reference_source_file(self, file_name):
        self.reference_star_list = []
        logging.info('Reading in table of reference sources...')
        table = np.loadtxt(file_name, skiprows=1)
        for row in table:
            x, y, flux = row
            self.reference_star_list.append(RefSource(x=y, y=x, f=flux))


    def apply_defaults(self):
        if 'file_list' not in self.__dict__:
            self.file_list = glob.glob(self.inDir + '*.fits')
        if 'star_list' not in self.__dict__:
            pass
        if 'reference_star_list' not in self.__dict__:
            if 'reference_star_file' not in self.__dict__:
                raise ValueError('Either a reference_star_file or a reference_star_list list must be provided to the PsfEstimator, but neither is given!')
            else:
                self.read_reference_source_file(self.reference_star_file)
        if 'tmpDir' not in self.__dict__:
            self.tmpDir = self.inDir + '../tmp/'
        if 'current_file' not in self.__dict__:
            setattr(self, 'current_file', 0)
        if 'current_frame' not in self.__dict__:
            setattr(self, 'current_frame', 0)
        if 'box_size' not in self.__dict__:
            setattr(self, 'box_size', 40)
        else:
            if isinstance(self.box_size, str):
                self.box_size = eval(self.box_size)
        if 'stats_box' not in self.__dict__:
            self.stats_box = Rectangle(None, None, None, None)
        else:
            print(type(eval(self.stats_box)))
            self.stats_box = Rectangle(**eval(self.stats_box))
            print(self.stats_box)
        if isinstance(self.noise_threshold, str):
            self.noise_threshold = eval(self.noise_threshold)


    def create_directories(self):
        for dir in ['tmpDir', 'outDir']:
            logging.info('Creating {} directory...'.format(dir))
            if not os.path.isdir(getattr(self, dir)):
                os.makedirs(getattr(self, dir))


    def initialize_psf_files(self):
        logging.info('Initializing psf files...')
        for file in self.file_list:
            logging.debug('Opening file {}...'.format(file))
            header = fits.getheader(file)
            psf_cube = np.zeros((header['NAXIS3'], self.box_size, self.box_size))
            psf_file = file.replace(self.inDir[:-1], self.tmpDir[:-1])
            psf_file = psf_file.replace('.fits', '_psf.fits')
            logging.info(file)
            logging.info(psf_file)
            fits.writeto(psf_file, data=psf_cube, overwrite=True)
        self.psf_file_list = glob.glob(self.tmpDir + '*_psf.fits')


    def estimate_psf(self, **kwargs):
        for file_index, file_name in enumerate(self.file_list):
            current_header = fits.getheader(file_name)
            number_of_frames = current_header['NAXIS3']
            for frame_index in range(number_of_frames):
                print('Estimating the PSF in file {:2}, frame {:4}.'.format(file_index + 1, frame_index + 1), end='\r')
                self.current_frame = fits.getdata(file_name)[frame_index]
                current_psf = self._estimate_psf(**kwargs)

                # save current psf to current psf file
                psf_file_name = self.psf_file_list[file_index]
                with fits.open(psf_file_name, mode='update') as hdulist:
                    hdulist[0].data[frame_index] = current_psf
                hdulist.flush()

        logging.info('Done: Estimating PSFs!')


    def _estimate_psf(self, **kwargs):
        if len(kwargs) is 0:

            background, noise = self.get_stats()

            mask = np.zeros((self.box_size, self.box_size))
            tmp = np.ma.masked_array(mask, mask=mask)
            for reference in self.reference_star_list:
                box = np.ma.masked_less(self.get_box(reference), background + self.noise_threshold * noise)
                # flux normalization
                box /= reference.flux
                # sum result of all reference stars to tmp
                tmp += box
            # normalization to unity
            tmp = tmp / np.sum(tmp)
            return tmp
        else:
            print('-----------------------s')

    def _estimate_psf_fine(self):
        """
        Estimate the psfs with a subgrid method, as used in the SUPERSTAR
        program by A. Marrasco et al.
        """
        print()



    def get_box(self, reference):
        radius = self.box_size // 2
        return self.current_frame[reference.x - radius : reference.x + radius + 1, reference.y - radius : reference.y + radius + 1]


    def get_stats(self):
        stats_box = self.stats_box(self.current_frame)
        return np.median(stats_box), np.std(stats_box)


    def reconstruction(self):
        """
        This function reconstructs the object following Eq. (1) from Schödel et al. (2013):
        O = {I_m P_m*} / {P_m P_m*}
        """
        pass
