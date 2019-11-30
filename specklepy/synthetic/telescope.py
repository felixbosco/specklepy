import os
import sys
import warnings
import numpy as np
from scipy.signal import fftconvolve
from scipy.ndimage import zoom
import astropy.units as u
from astropy.io import fits
from astropy.modeling import models



class Telescope(object):

	"""Class carrying the information of a telescope.

	Attributes:
		diameter (astropy.units.Quantity):
		psf_source (str):
		psf_plane (int):

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


	def __init__(self, diameter, psf_source, central_obscuration=None, psf_plane=0, **kwargs):
		"""Instantiate Telescope class:

		Args:
			diameter (astrop.units.Quantity): Telescope diameter, used to
				compute the light collecting area.
			psf_source (str): File name to read PSFs from or model name. Models
				can be either 'AiryDisk' or 'Gaussian'. The models require
			central_obscuration (float, optional): Radial fraction of the
				telescope aperture that is blocked by the secondary.
			psf_plane (int, optional): Index of the first frame to read from
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

		if isinstance(psf_plane, int):
			self.psf_plane = psf_plane
		else:
			raise TypeError(self.typeerror.format('psf_plane', type(psf_plane), 'int'))

		# Derive secondary parameters
		if self.central_obscuration is not None:
			self.area = (1. - self.central_obscuration**2) * np.pi * (self.diameter / 2)**2
		else:
			self.area = np.pi * (self.diameter / 2)**2

		if psf_source.lower() in ['airydisk', 'gaussian']:
			self.model_psf(psf_source, **kwargs)
		else:
			self.read_psf_file(psf_source)



	def __call__(self, flux_array, flux_array_resolution, integration_time=None, verbose=0, **kwargs):
		"""Requires the integration_time only if the PSF is non-static."""
		tmp = flux_array * self.area
		total_flux = np.sum(tmp)
		tmp_unit = tmp.unit

		# Resample flux_array to psf resolution
		try:
			ratio = float(flux_array_resolution / self.psf_resolution)
			#if ratio < 1.0:
			#	print('')
			#	print('Be cautious, the resolution of the PSF is worse than of the target data.')
			#	#raise ValueError('Be cautious, the resolution of the PSF  ({}) is worse than of the target data ({}).'.format(self.psf_resolutio, flux_array_resolution))
		except UnitConversionError as e:
			raise UnitConversionError("The resolution values of the image ({}) and PSF ({}) have different units!".format(flux_array_resolution, self.psf_resolution))

		# Prepare PSF if non-static
		if hasattr(self, 'timestep'):
			if integration_time is None:
				raise ValueError("If the PSF source of Telescope is non-static, the call function requires the integration_time.")
			self.integrate_psf(integration_time=integration_time)

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
		if tmp.shape[0] > 2048+512 or self.psf.shape[0] > 512+256:
			print('With these sizes (image: {} and PSF: {}), the computation will be very expensive. It is suggested to adapt the resolution of the objects.'.format(tmp.shape, self.psf.shape))
			user_input = input('Do you still want to continue? [Y/N]')
			if user_input in ['Y', 'y', 'Yes', 'yes']:
				pass
			else:
				raise Exception('Program aborted, re-define the resolution of the objects.')
		convolved = fftconvolve(tmp, self.psf, mode='same') * tmp_unit
		if verbose > 0:
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



	def _get_value(self, header, key, alias_dict={'sec': 's', 'milliarcsec': 'mas', 'microns': 'micron'}, verbose=False):
		"""
		The alias_dict dictionary is used as a mapping from
		unvalid unit strings for u.Unit(str). Feel free to
		add new aliases.
		"""

		value = header[key]
		unit_str = header.comments[key]

		if unit_str == '':
		    if verbose:
		        print("Function 'get_value()' did not find a unit in the comment string.")
		    return value
		else:
		    try:
		        unit = u.Unit(unit_str)
		    except ValueError as e:
		        if verbose:
		            print("ValueError:", e)
		            print("Trying aliases from alias_dict...")
		        try:
		            unit = u.Unit(alias_dict[unit_str])
		        except:
		            raise IOError("Found no matching key in the alias_dict. You may add the corresponding entry.")
		    return value * unit



	def normalize(self, array, mode='unity_circular'):
		"""Normalizes the array to either have a sum of 1 ('unity' mode) or that the peak value is 1 ('max' mode)."""
		if mode == 'unity':
		    return array / np.sum(array)
		elif mode == 'max':
		    return array / np.max(array)
		elif mode == 'unity_circular':
		    x, y = array.shape
		    low_cut = array[0, int(y/2)]
		    array = np.maximum(array - low_cut, 0)
		    return array / np.sum(array)



	def integrate_psf(self, integration_time, hdu_entry=0):
		number_planes = int(integration_time / self.timestep)
		with fits.open(self.psf_source) as hdulist:
			data = hdulist[hdu_entry].data

			self.psf_plane += 1
			if self.psf_plane + number_planes < data.shape[0]:
				self.psf = np.sum(data[self.psf_plane : self.psf_plane+number_planes], axis=0)
			else:
				self.psf = np.sum(data[self.psf_plane : ], axis=0)
				self.psf += np.sum(data[ : (self.psf_plane+number_planes) % data.shape[0]], axis=0)
			self.psf_plane += number_planes - 1
			self.psf_plane = self.psf_plane % data.shape[0]
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
			model = models.AiryDisk2D(x_0=center[0], y_0=center[1], radius=float(radius / psf_resolution))
		elif model.lower() == 'gaussian':
			model = models.Gaussian2D(x_mean=center[0], y_mean=center[1], x_stddev=float(radius / psf_resolution), y_stddev=float(radius / psf_resolution))
		else:
			raise ValueError("model_psf received model argument {}, but must be\
							either 'AriyDisk' or 'Gaussian'!".format(model))

		y, x = np.mgrid[0:shape[0], 0:shape[1]]
		self.psf = model(x, y)
		self.psf = self.normalize(self.psf)
