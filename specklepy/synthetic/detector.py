import sys
import warnings
import numpy as np
from scipy.ndimage import zoom
from astropy import units as u

from specklepy.logging import logging
from specklepy.utils.plot import imshow



class Detector(object):

	"""Class carrying information of an astronomical detector.

	Attributes:
		shape (tuple, dtype=int): 2D tuple of the number of pixels.
		pixel_scale (astropy.units.Quantity): Pixel scale in arcsec.

	Optional attributes are:
		quantum_efficiency (u.Quantity):
		system_gain (u.Quantity):
		readout_noise (u.Quantity):
		dark_current (u.Quantity):
		saturation_level (u.Quantity):
		optics_transmission (u.Quantity):

	Future features:
		Attribute dictionary 'readout_modes': Shall enable the flexible use of different
			readout modes with different parameters.
		Pixel non-uniformity: Shall use a flatfield file to create static
			non-uniform pixel responses.
	"""

	__name__ = 'detector'
	typeerror = 'Detector received {} argument of {} type, but needs to be {}!'



	def __init__(self, shape, pixel_scale, optics_transmission=1, quantum_efficiency=1, system_gain=1, readout_noise=0, dark_current=None, saturation_level=None):
		"""Instantiate Detector class.

		Args:
			shape (tuple, dtype=int): Shape of the detector array, i.e. number
				of pixels. If provided as int, then the detector will be square
				shaped.
			pixel_scale (u.Quantity):
			optics_transmission (u.Quantity, optional):
			quantum_efficiency (u.Quantity, optional):
			system_gain (u.Quantity, optional):
			readout_noise (u.Quantity, optional):
			dark_current (u.Quantity, optional):
			saturation_level (u.Quantity, optional):
		"""

		# Input parameters
		if isinstance(shape, tuple):
			self.shape = shape
		elif isinstance(shape, int):
			self.shape = (shape, shape)
		else:
			raise TypeError(self.typeerror.format('shape', type(shape), 'tuple'))

		if isinstance(pixel_scale, u.Quantity):
			self.pixel_scale = pixel_scale
		elif isinstance(pixel_scale, float) or isinstance(pixel_scale, int):
			logging.warning("Interpreting float type pixel_scale as {}".format(pixel_scale * u.arcsec))
			self.pixel_scale = pixel_scale * u.arcsec
		else:
			raise TypeError(self.typeerror.format('pixel_scale', type(pixel_scale), 'u.Quantity'))

		if isinstance(quantum_efficiency, u.Quantity):
			self.quantum_efficiency = quantum_efficiency
		elif isinstance(quantum_efficiency, float) or isinstance(quantum_efficiency, int):
			logging.warning("Interpreting float type quantum_efficiency as {}".format(quantum_efficiency * u.electron / u.ph))
			self.quantum_efficiency = quantum_efficiency * u.electron / u.ph
		else:
			raise TypeError(self.typeerror.format('quantum_efficiency', type(quantum_efficiency), 'u.Quantity'))

		if isinstance(system_gain, u.Quantity):
			self.system_gain = system_gain
		elif isinstance(system_gain, float) or isinstance(system_gain, int):
			logging.warning("Interpreting float type system_gain as {}".format(system_gain * u.electron / u.adu))
			self.system_gain = system_gain * u.electron / u.adu
		else:
			raise TypeError(self.typeerror.format('system_gain', type(system_gain), 'u.Quantity'))

		if isinstance(readout_noise, u.Quantity):
			self.readout_noise = readout_noise
		elif isinstance(readout_noise, float) or isinstance(readout_noise, int):
			logging.warning("Interpreting float type readout_noise as {}".format(readout_noise * u.electron))
			self.readout_noise = readout_noise * u.electron
		else:
			raise TypeError(self.typeerror.format('readout_noise', type(readout_noise), 'u.Quantity'))

		if isinstance(optics_transmission, float) or isinstance(optics_transmission, int):
			self.optics_transmission = optics_transmission
		else:
			raise TypeError(self.typeerror.format('optics_transmission', type(optics_transmission), 'float'))

		if isinstance(saturation_level, u.Quantity) or saturation_level is None:
			self.saturation_level = saturation_level
		elif isinstance(saturation_level, float) or isinstance(saturation_level, int):
			logging.warning("Interpreting float type saturation_level as {}".format(saturation_level * u.electron))
			self.saturation_level = saturation_level * u.electron
		else:
			raise TypeError(self.typeerror.format('saturation_level', type(saturation_level), 'u.Quantity'))

		# Derive secondary parameters
		self.array = np.zeros(self.shape)
		self.FoV = (self.shape[0] * self.pixel_scale, self.shape[1] * self.pixel_scale)



	@property
	def resolution(self):
		return self.pixel_scale

	@resolution.setter
	def resolution(self, value):
		self.pixel_scale = value



	def __call__(self, *args, **kwargs):
		return self.get_counts(*args, **kwargs)



	def __str__(self):
		tmp = "Detector:\n"
		for key in self.__dict__:
			if key == 'array':
				continue
			tmp += "{}: {}\n".format(key, self.__dict__[key])
		return tmp



	def get_counts(self, photon_rate, integration_time, target_FoV, compute_photon_shot_noise=True, debug=False):
		self.expose(photon_rate, integration_time, target_FoV, compute_photon_shot_noise=compute_photon_shot_noise, debug=debug)
		return self.readout(integration_time=integration_time)



	def resample(self, target_data, target_FoV):
		if target_FoV[0] < self.FoV[0] or target_FoV[1] < self.FoV[1]:
			raise ValueError('The FoV of the target object ({}) is smaller than that of the detector ({})!'.format(target_FoV, self.FoV))
		subfield_shape = [round((self.FoV[i]/ target_FoV[i]).decompose().value * target_data.shape[i]) for i, val in enumerate(target_FoV)]
		x0 = int( (target_data.shape[0] - subfield_shape[0])/2 )
		y0 = int( (target_data.shape[1] - subfield_shape[1])/2 )
		subfield = target_data[x0:x0+subfield_shape[0], y0:y0+subfield_shape[1]]
		stretch_ratio = (self.shape[0] / subfield.shape[0], self.shape[1] / subfield.shape[1])
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			return zoom(subfield, stretch_ratio, order=1) / stretch_ratio[0] / stretch_ratio[1] * subfield.unit



	def expose(self, photon_rate, integration_time, target_FoV, compute_photon_shot_noise=True, debug=False):
		if debug:
			imshow(photon_rate, title='photon_rate')
		tmp = self.resample(photon_rate, target_FoV) * integration_time
		tmp *= self.quantum_efficiency
		if hasattr(self, 'optics_transmission'):
			tmp *= self.optics_transmission
		if compute_photon_shot_noise:
			if debug:
				imshow(tmp, title='expose : tmp')
			try:
				tmp = np.random.poisson(tmp.value) * tmp.unit
			except ValueError as e:
				if debug:
					print('Bypassed ValueError ({}) in np.random.poisson() by substituting values smaller than zero by zero.'.format(e))
				tmp = np.random.poisson(np.maximum(tmp.value, 0.0)) * tmp.unit
		if self.saturation_level is not None:
			tmp = np.minimum(tmp, self.saturation_level) # * self.system_gain)
		self.array = np.round(tmp)



	def readout(self, integration_time, reset=True, window=None):
        # Read copy and clear the array
		tmp = self.array
		if hasattr(self, 'dark_current'):
			tmp += np.random.poisson(self.dark_current.value, self.shape) * self.dark_current.unit * u.pix * integration_time
		if hasattr(self, 'readout_noise'):
			tmp += np.round(np.random.normal(0.0, self.readout_noise.value, self.shape) ) * self.readout_noise.unit * u.pix
		tmp /= self.system_gain
		if self.saturation_level is not None:
			return np.minimum(tmp, self.saturation_level / self.system_gain)
		# Reset the detector
		if reset:
			self.array = np.zeros(self.shape)
		return tmp.decompose()
