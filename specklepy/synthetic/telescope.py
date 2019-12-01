import os
import sys
import warnings
import numpy as np
from scipy.signal import fftconvolve
from scipy.ndimage import zoom
import astropy.units as u
from astropy.io import fits
from astropy.modeling import models

from specklepy.utils.plot import imshow
from specklepy.logging import logging



class Telescope(object):

	"""Class carrying the information of a telescope.

	Attributes:
		diameter (astropy.units.Quantity):
		psf_source (str):
		psf_frame (int):

	Optional attributes:
		central_obscuration (float, optional):

	Future features:
		Mode 'airy_model': Shall compute the psf with the Fourier transform of
			the aperture instead of a file.
		Mode 'seeing': Shall compute the psf as a Gaussian seeing disk instead
			of a file.
	"""

	__name__ = 'telescope'
	typeerror = 'Telescope received {} argument of {} type, but needs to be {}!'
	TIME_STEP_KEYS = ['TIMESTEP', 'INTTIME', 'CDELT3']
	RESOLUTION_KEYS = ['PIXSIZE', 'CDELT1']


	def __init__(self, diameter, psf_source, central_obscuration=None, psf_frame=0, **kwargs):
		"""Instantiate Telescope class:

		Args:
			diameter (astrop.units.Quantity): Telescope diameter, used to
				compute the light collecting area.
			psf_source (str): File name to read PSFs from or model name. Models
				can be either 'AiryDisk' or 'Gaussian'. The models require
			central_obscuration (float, optional): Radial fraction of the
				telescope aperture that is blocked by the secondary.
			psf_frame (int, optional): Index of the first frame to read from
				psf_source.
			kwargs: Are forwarded to the psf_source model.
		"""

		# Input parameters
		if isinstance(diameter, u.Quantity):
			self.diameter = diameter
		elif isinstance(diameter, float) or isinstance(diameter, int):
			logging.warning("Interpreting float type diameter as {}".format(diameter * u.m))
			self.diameter = diameter * u.m
		else:
			raise TypeError(self.typeerror.format('diameter', type(diameter), 'u.Quantity'))

		if isinstance(psf_source, str):
			self.psf_source = psf_source
		else:
			raise TypeError(self.typeerror.format('psf_source', type(psf_source), 'str'))

		if isinstance(central_obscuration, float) or central_obscuration is None:
			self.central_obscuration = central_obscuration
		else:
			raise TypeError(self.typeerror.format('central_obscuration', type(central_obscuration), 'float'))

		if isinstance(psf_frame, int):
			self.psf_frame = psf_frame
		else:
			raise TypeError(self.typeerror.format('psf_frame', type(psf_frame), 'int'))

		# Derive secondary parameters
		if self.central_obscuration is not None:
			self.area = (1. - self.central_obscuration**2) * np.pi * (self.diameter / 2)**2
		else:
			self.area = np.pi * (self.diameter / 2)**2

		if psf_source.lower() in ['airydisk', 'gaussian']:
			self.model_psf(psf_source, **kwargs)
		else:
			self.read_psf_file(psf_source)


	def __call__(self, *args, **kwargs):
		return self.get_photon_rate(*args, **kwargs)


	def get_photon_rate(self, photon_rate_density, photon_rate_density_resolution=None, integration_time=None, debug=False):
		"""Propagates the 'photon_rate_density' array through the telescope.

		The photon_rate_density is multiplied by the telescope collecting area and then
		convolved with the PSF. If the resolution of the flux array is different
		from the telescopes psf_resolution, then one is resampled. If the PSF is
		non-static, it will be integrated over the 'integration_time' value.

		Args:
			photon_rate_density (np.ndarray, dtype=u.Quantity):
			photon_rate_density_resolution (u.Quantity, optional):
			integration_time(u.Quantity, optional): Required only if the PSF is
				non-static.
			debug (bool, optional): Set True for debugging. Default is False.

		Returns:
			photon_rate (u.Quantity): PSF-convolved photon rate array.
		"""

		# Input parameters
		if not isinstance(photon_rate_density, u.Quantity):
			raise TypeError(self.typeerror.format('photon_rate_density', type(photon_rate_density), 'u.Quantity'))

		if photon_rate_density_resolution is not None:
			if not isinstance(photon_rate_density_resolution, u.Quantity):
				raise TypeError(self.typeerror.format('photon_rate_density_resolution', type(photon_rate_density_resolution), 'u.Quantity'))
			psf_resample_mode = True
		else:
			psf_resample_mode = False

		if integration_time is None and hasattr(self, 'timestep'):
			raise ValueError("If the PSF source of Telescope is non-static, the call function requires the integration_time.")
		elif isinstance(integration_time, float) or isinstance(integration_time, int):
			logging.warning("Interpreting float type integration_time as {}".format(integration_time * u.s))
			integration_time = integration_time * u.s
		elif not isinstance(integration_time, u.Quantity):
			raise TypeError(self.typeerror.format('integration_time', type(integration_time), 'u.Quantity'))


		# Apply telescope collecting area
		tmp = photon_rate_density * self.area
		total_flux = np.sum(tmp)
		tmp_unit = tmp.unit


		# Prepare PSF if non-static
		if hasattr(self, 'timestep'):
			self.integrate_psf(integration_time=integration_time)


		# Resample photon_rate_density to psf resolution
		if psf_resample_mode:
			try:
				ratio = float(photon_rate_density_resolution / self.psf_resolution)
			except UnitConversionError as e:
				raise UnitConversionError("The resolution values of the image ({}) and PSF ({}) have different units!".format(photon_rate_density_resolution, self.psf_resolution))

			#convolved = np.zeros(tmp.shape)
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				if ratio < 1.0:
					self.psf = zoom(self.psf, 1/ratio, order=1) / ratio**2
					self.psf = self.normalize(self.psf)
				else:
					memory_sum = np.sum(tmp)
					tmp = zoom(tmp, ratio, order=1) / ratio**2
					tmp = tmp / np.sum(tmp) * memory_sum
			# if tmp.shape[0] > 2048+512 or self.psf.shape[0] > 512+256:
			# 	print('With these sizes (image: {} and PSF: {}), the computation will be very expensive. It is suggested to adapt the resolution of the objects.'.format(tmp.shape, self.psf.shape))
			# 	user_input = input('Do you still want to continue? [Y/N]')
			# 	if user_input in ['Y', 'y', 'Yes', 'yes']:
			# 		pass
			# 	else:
			# 		raise Exception('Program aborted, re-define the resolution of the objects.')
		convolved = fftconvolve(tmp, self.psf, mode='same') * tmp_unit
		if debug:
			print('Check of flux conservation during convolution:')
			print('Before: ', total_flux)
			print('After:  ', np.sum(convolved))
		return convolved.decompose()



	def __str__(self):
		tmp = "Telescope:\n"
		for key in self.__dict__:
			if key == 'psf':
				continue
			tmp += "{}: {}\n".format(key, self.__dict__[key])
		return tmp



	def read_psf_file(self, filename, hdu_entry=0):
		with fits.open(filename) as hdulist:
			header = hdulist[hdu_entry].header
		if header['NAXIS'] == 2:
			with fits.open(self.psf_source) as hdulist:
				self.psf = self.normalize(hdulist[hdu_entry].data)
		else:
			for key in self.TIME_STEP_KEYS:
				try:
					self.timestep = self._get_value(header, key)
					break
				except KeyError as e:
					continue
		for key in self.RESOLUTION_KEYS:
			try:
				self.psf_resolution = self._get_value(header, key)
				break
			except KeyError as e:
				continue
			raise IOError("No key from {} was found in file for the psf resolution.".format(self.RESOLUTION_KEYS))



	def _get_value(self, header, key, alias_dict={'sec': 's', 'milliarcsec': 'mas', 'microns': 'micron'}, debug=False):
		"""
		The alias_dict dictionary is used as a mapping from
		unvalid unit strings for u.Unit(str). Feel free to
		add new aliases.
		"""

		value = header[key]
		unit_str = header.comments[key]

		if unit_str == '':
		    if debug:
		        print("Function 'get_value()' did not find a unit in the comment string.")
		    return value
		else:
		    try:
		        unit = u.Unit(unit_str)
		    except ValueError as e:
		        if debug:
		            print("ValueError:", e)
		            print("Trying aliases from alias_dict...")
		        try:
		            unit = u.Unit(alias_dict[unit_str])
		        except:
		            raise IOError("Found no matching key in the alias_dict. You may add the corresponding entry.")
		    return value * unit



	def normalize(self, array, mode='sum_circular'):
		"""Normalizes the input array depending on the mode.

		Args:
			array (np.ndarray): Array to be normalized.
			mode (str, optional): Can be either 'sum' for having a sum of 1,
				'peak' for having a peak value 1, or 'sum_circular' for
				subtracting a constant and then normalizing to a sum of 1.
				Default is 'sum_circular'.

		Returns:
			Normalized array (np.ndarray)
		"""

		if not isinstance(array, np.ndarray):
			raise TypeError(self.typeerror.format('array', type(array), 'np.ndarray'))
		if np.sum(array) == 0:
			raise ValueError("Normalize received an array of zeros!")

		if mode == 'sum':
		    return array / np.sum(array)
		elif mode == 'max':
		    return array / np.max(array)
		elif mode == 'sum_circular':
		    x, y = array.shape
		    low_cut = array[0, int(y/2)]
		    array = np.maximum(array - low_cut, 0)
		    return self.normalize(array, mode='sum')



	def integrate_psf(self, integration_time, hdu_entry=0, debug=False):
		"""Integrates psf frames over the input time.

		Args:
			integration_time (u.Quantity): This is used to compute number of
				frames 'nframes', via floor division by the timestep attribute.
				If None, then the function just exits.
			debug (bool, optional): Set True for debugging. Default is False.
		"""

		if isinstance(integration_time, int) or isinstance(integration_time, float):
			logging.warning("Interpreting float type integration_time as {}".format(integration_time * u.s))
			integration_time = integration_time * u.s
		elif not isinstance(integration_time, u.Quantity):
			raise TypeError('integrate_psf received integration_time argument of type {}, but needs to be u.Quantity')

		if integration_time < self.timestep:
			raise ValueError("The integration time {} was chosen shorter than the time resolution of the psf source, of {}".format(integration_time, self.timestep))

		nframes = int(integration_time / self.timestep)

		with fits.open(self.psf_source) as hdulist:
			data = hdulist[hdu_entry].data

			self.psf_frame += 1
			if self.psf_frame + nframes < data.shape[0]:
				self.psf = np.sum(data[self.psf_frame : self.psf_frame+nframes], axis=0)
			else:
				self.psf = np.sum(data[self.psf_frame : ], axis=0)
				self.psf += np.sum(data[ : (self.psf_frame+nframes) % data.shape[0]], axis=0)
			self.psf_frame += nframes - 1
			self.psf_frame = self.psf_frame % data.shape[0]

			#Normalization
			self.psf = self.normalize(self.psf)



	def model_psf(self, model, radius, psf_resolution, shape=256, **kwargs):
		"""Models the PSF given the desired model function and kwargs.

		Args:
			model (str): Must be either 'airydisk' or 'gaussian'.
			kwargs are forwarded to the model function.
		"""

		if not isinstance(model, str):
			raise TypeError('model_psf received model argument of type {}, but \
							needs to be str type!')

		if isinstance(radius, u.Quantity):
			self.radius = radius
		elif isinstance(radius, float) or isinstance(radius, int):
			logging.warning("Interpreting float type radius as {}".format(radius * u.arcsec))
			self.radius = radius * u.arcsec
		else:
			raise TypeError(self.typeerror.format('radius', type(radius), 'u.Quantity'))

		if isinstance(psf_resolution, u.Quantity):
			self.psf_resolution = psf_resolution
		elif isinstance(psf_resolution, float) or isinstance(psf_resolution, int):
			logging.warning("Interpreting float type psf_resolution as {}".format(psf_resolution * u.arcsec))
			self.psf_resolution = psf_resolution * u.arcsec
		else:
			raise TypeError(self.typeerror.format('psf_resolution', type(psf_resolution), 'u.Quantity'))

		if isinstance(shape, int):
			center = (shape / 2, shape / 2)
			shape = (shape, shape)
		elif isinstance(shape, tuple):
			center = (shape[0] / 2, shape[1] / 2)
		else:
			raise TypeError('model_psf received shape argument of type {}, but \
							needs to be int or tuple type!')

		if model.lower() == 'airydisk':
			model = models.AiryDisk2D(x_0=center[0], y_0=center[1], radius=float(self.radius / self.psf_resolution))
		elif model.lower() == 'gaussian':
			model = models.Gaussian2D(x_mean=center[0], y_mean=center[1], x_stddev=float(self.radius / self.psf_resolution), y_stddev=float(self.radius / self.psf_resolution))
		else:
			raise ValueError("model_psf received model argument {}, but must be\
							either 'AriyDisk' or 'Gaussian'!".format(model))

		y, x = np.mgrid[0:shape[0], 0:shape[1]]
		self.psf = model(x, y)
		self.psf = self.normalize(self.psf)
