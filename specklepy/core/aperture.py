from copy import copy
import numpy as np
import warnings

from specklepy.io import fits
from specklepy.logging import logger
from specklepy.utils.box import Box, Box3D
from specklepy.utils import transferfunctions as tf


class Aperture(object):

    variance_extension = 'VAR'
    mask_extension = 'MASK'

    def __init__(self, *args, file_name, path=None, mask='circular', crop=True, fill_value=0):
        """Instantiate Aperture object.

        Long description...

        Args:
            *args (int, float):
                Can be either two or three arguments, where the last is always interpreted as the aperture radius,
                which must be int type. The other arguments can be either a coordinate tuple or two individual
                coordinates.
            file_name (np.ndarray, str):
                2D or 3D np.ndarray that the aperture shall be extracted of. If provided as str type, this is assumed
                to be a path and the objects tries to read the FITS file.
            mask (str, optional):
                Mode that is describing, how the aperture is masked. Can be 'circular' or 'rectangular'. If 'circular',
                then it creates a circular mask and the data become a np.ma.masked_array. Default is 'circular'.
            crop (bool, optional):
                If set to True, then the object only stores a copy of the data with radius around the center. Otherwise
                all the data beyond the limits of the aperture are masked. Default is True.
            fill_value (int or float, optional):
                Value for filling masked pixels. Default is `0`.
        """

        # Interpret the args
        if len(args) == 2 and (isinstance(args[0], tuple) or isinstance(args[0], list)):
            y0 = args[0][0]
            x0 = args[0][1]
            radius = args[1]
        elif len(args) == 3:
            y0 = args[0]
            x0 = args[1]
            radius = args[2]
        else:
            raise ValueError(f"Aperture expects either 2 or 3 positional args, but {len(args)} have been provided!")

        # Integer checks and handling non-integer input
        if isinstance(y0, int) and isinstance(x0, int):
            self.y0 = y0
            self.x0 = x0
            self.y_offset = 0
            self.x_offset = 0
        else:
            y0 = float(y0)
            x0 = float(x0)
            self.y0 = np.rint(y0).astype(int)
            self.x0 = np.rint(x0).astype(int)
            self.y_offset = y0 - self.y0
            self.x_offset = x0 - self.x0

        if isinstance(radius, int):
            self.radius = radius
        else:
            raise ValueError(f"Aperture radius must be integer type, but is {radius})")

        self.fill_value = fill_value

        # Handling data input
        if isinstance(file_name, str):
            logger.debug(f"Aperture argument data '{file_name}' is interpreted as file name.")
            self.data = fits.get_data(file_name)
            self.var = fits.get_data(file_name=file_name, path=path, extension=self.variance_extension,
                                     ignore_missing_extension=True)
            data_mask = fits.get_data(file_name=file_name , path=path, extension=self.mask_extension, dtype=bool,
                                      ignore_missing_extension=True)
            if data_mask is not None:
                self.data[data_mask] = self.fill_value

        else:
            raise TypeError

        # Remove the margins if requested
        self.cropped = False
        if crop:
            self.crop()

        # Create a mask
        self.mask_mode = mask
        mask = self.make_mask(mode=self.mask_mode)
        self.data = np.ma.masked_array(self.data, mask=mask)

        # Initialize optional attributes
        self.power_spectrum_cube = None

    @property
    def width(self):
        return self.radius * 2 + 1

    @property
    def index(self):
        return self.y0, self.x0

    @property
    def offset(self):
        return self.y_offset, self.x_offset

    @property
    def shape(self):
        return self.data.shape

    @property
    def n_frames(self):
        return self.data.shape[0]

    def __getitem__(self, index):
        return self.data[index]

    def make_radius_map(self):
        if self.cropped:
            center = (self.radius, self.radius)
        else:
            center = self.index
        xx, yy = np.mgrid[:self.data.shape[-2], :self.data.shape[-1]]
        return np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))

    def make_mask(self, mode='circular'):
        """Create a circular or rectangular mask."""

        if not isinstance(mode, str):
            raise TypeError(f"Aperture received mode argument of type {type(mode)}, but needs to be str type!")

        if mode == 'circular':
            distance_map = self.make_radius_map()
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
                    gpm = Box.centered_at(self.x0, y0=self.y0, radius=self.radius)
                    mask = gpm.set_value(mask, False)
                elif mask.ndim == 3:
                    gpm = Box3D.centered_at(self.x0, y0=self.y0, z0=None, radius=self.radius)
                    mask = gpm.set_value(mask, False)
                return mask
        else:
            raise ValueError(f"Aperture received mode argument {mode}, but needs to be either 'circular' or "
                             f"'rectangular'!")

    def crop(self):
        if self.cropped:
            logger.info("Margins are removed already from aperture instance.")
        else:
            if self.data.ndim == 2:
                gpm = Box.centered_at(self.x0, y0=self.y0, radius=self.radius)
                self.data = copy(gpm(self.data))
            elif self.data.ndim == 3:
                gpm = Box3D.centered_at(self.x0, y0=self.y0, z0=None, radius=self.radius)
                self.data = copy(gpm(self.data))
            self.cropped = True

    def remove_margins(self):
        self.crop()

    def get_integrated(self):
        """Returns a 2-dimensional representation of the aperture and integrates
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

    def center_on_peak(self):
        # Make sure that the margins around the aperture are available
        if self.cropped:
            raise NotImplementedError("Centering cropped apertures is not implemented yet!")

        # Estimate intensity peak
        peak = self.get_aperture_peak()
        logger.info(f"Re-centering the aperture from {self.index} to {peak}")

        # Update properties
        self.y0, self.x0 = peak
        self.y_offset, self.x_offset = (0, 0)

        # Update mask
        mask = self.make_mask(mode=self.mask_mode)
        self.data = np.ma.masked_array(self.data, mask=mask)

    def get_power_spectrum_cube(self):
        self.power_spectrum_cube = np.zeros(self.data.shape)
        for index, frame in enumerate(self.data):
            self.power_spectrum_cube[index] = tf.powerspec(frame)
        logger.info("Computed the Fourier transform of every frame.")

    def get_power_spectrum_profile(self):

        if self.power_spectrum_cube is None:
            self.get_power_spectrum_cube()

        # Initialize output radii and array
        radius_map = self.make_radius_map()
        rdata = np.unique(radius_map)
        ydata = np.zeros(rdata.shape)
        edata = np.zeros(rdata.shape)

        # Iterate over aperture radii
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for index, radius in enumerate(rdata):
                y, x = np.where(radius_map == radius)
                subset = self.power_spectrum_cube[:, y, x]
                ydata[index] = np.mean(subset)
                edata[index] = np.std(subset)

        return rdata, ydata, edata

    def get_psf_variance(self):
        """Extract a radial variance profile of the speckle PSFs."""

        # Check data shape
        if self.data.ndim != 3:
            logger.warning(f"Aperture data need to be 3D but are {self.data.ndim}-dimensional!")
            raise SystemExit("Aborting...")

        # Initialize output radii and array
        radius_map = self.make_radius_map()
        rdata = np.unique(radius_map)
        ydata = np.zeros(rdata.shape)
        edata = np.zeros(rdata.shape)

        # Extract 2D variance map
        var_map = np.var(self.data, axis=0)

        # Iterate over aperture radii
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for index, radius in enumerate(rdata):
                subset = var_map[np.where(radius_map == radius)]
                ydata[index] = np.sqrt(np.mean(subset))
                edata[index] = np.sqrt(np.var(subset))

        return rdata, ydata, edata

    def get_encircled_energy(self, saveto=None):
        """Extracts the encircled energy from an aperture as a function of
        radius."""

        # Initialize output radii and array
        radius_map = self.make_radius_map()
        rdata = np.unique(radius_map)
        ydata = np.zeros(rdata.shape)
        edata = np.zeros(rdata.shape)

        # Extract 2D image
        image = self.get_integrated()

        # Iterate over aperture radii
        for index, radius in enumerate(rdata):
            subset = image[np.where(radius_map <= radius)]
            ydata[index] = np.sum(subset)
            edata[index] = np.std(subset)

        # Save results to file
        if saveto is not None:
            header = "radius_(pix) encircled_energy(data_unit)"
            data = np.concatenate(([rdata], [ydata]), axis=0).transpose()
            np.savetxt(saveto, data, header=header)
            logger.info(f"Saved encircled energy data to file {saveto}")

        return rdata, ydata, edata

    def get_psf_profile(self):
        """Computes the radial average of the PSF in the aperture."""

        # Initialize output radii and array
        radius_map = self.make_radius_map()
        rdata = np.unique(radius_map)
        ydata = np.zeros(rdata.shape)
        edata = np.zeros(rdata.shape)

        # Extract 2D image
        image = self.get_integrated()

        # Iterate over aperture radii
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for index, radius in enumerate(rdata):
                subset = image[np.where(radius_map == radius)]
                ydata[index] = np.mean(subset)
                edata[index] = np.std(subset)

        return rdata, ydata, edata

    def spatial_frequency(self, pixel_scale=1, profile=True):
        r = self.make_radius_map()
        if profile:
            return np.divide(np.unique(r), self.radius) / (2 * pixel_scale)
        else:
            return np.divide(r, self.radius) / (2 * pixel_scale)

    def spatial_wavelength(self, pixel_scale=1, profile=True):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.divide(1, self.spatial_frequency(pixel_scale=pixel_scale, profile=profile))
