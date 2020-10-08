import numpy as np
from scipy.ndimage import zoom
import warnings

from astropy.units import Unit, Quantity

from specklepy.exceptions import SpecklepyTypeError
from specklepy.io import config
from specklepy.logging import logger
from specklepy.utils.plot import imshow


class Detector(object):

	"""Class carrying information of an astronomical detector.

	Attributes:
		shape (tuple, dtype=int):
			Tuple of the number of pixels along two axes.
		pixel_scale (astropy.units.Quantity):
			Angular size of a pixel in arcsec.

	Optional attributes are:
		quantum_efficiency (Quantity):
		system_gain (Quantity):
		readout_noise (Quantity):
		dark_current (Quantity):
		saturation_level (Quantity):
		optics_transmission (Quantity):

	Future features:
		Attribute dictionary 'readout_modes':
			Shall enable the flexible use of different readout modes with different parameters.
		Pixel non-uniformity:
			Shall use a flatfield file to create static non-uniform pixel responses.
	"""

	__name__ = 'detector'

	def __init__(self, shape, pixel_scale, optics_transmission=1, quantum_efficiency=1, system_gain=1,
				 readout_noise=None, dark_current=None, saturation_level=None):
		"""Instantiate Detector class.

		Args:
			shape (tuple, dtype=int):
				Number of the pixels of the detector along two axes. Integer values will create a square detector array
				with the same number of pixels along both axes.
			pixel_scale (Quantity):
				Angular size of each pixel.
			optics_transmission (float, optional):
				Optical transmission coefficient, scaling between 0.0 and 1.0 (default).
			quantum_efficiency (float, Quantity, optional):
				Quantum efficiency of the detector. Scalar type values will be interpreted units of electrons per
				photon.
			system_gain (float, Quantity, optional):
				System gain of the detector. Scalar type values will be interpreted units of electrons per ADU.
			readout_noise (float, Quantity, optional):
				Read noise of the detector. Scalar type values will be interpreted units of electrons.
			dark_current (float, Quantity, optional):
				Dark current of the detector. Scalar type values will be interpreted units of electrons per second.
			saturation_level (float, Quantity, optional):
				Saturation level of the detector. Scalar type values will be interpreted units of electrons.
		"""

		# Input parameters
		if isinstance(shape, str):
			shape = eval(shape)
		if isinstance(shape, tuple):
			self.shape = shape
		elif isinstance(shape, int):
			self.shape = (shape, shape)
		else:
			raise SpecklepyTypeError('Detector', 'shape', type(shape), 'tuple')

		if isinstance(pixel_scale, Quantity):
			self.pixel_scale = pixel_scale
		elif isinstance(pixel_scale, (int, float)):
			logger.warning(f"Interpreting float type pixel_scale as {pixel_scale} arcsec")
			self.pixel_scale = pixel_scale * Unit('arcsec')
		elif isinstance(pixel_scale, str):
			self.pixel_scale = Quantity(pixel_scale)
		else:
			raise SpecklepyTypeError('Detector', 'pixel_scale', type(pixel_scale), 'Quantity')

		if isinstance(optics_transmission, (int, float)):
			self.optics_transmission = optics_transmission
		elif isinstance(optics_transmission, str):
			self.optics_transmission = float(optics_transmission)
		else:
			raise SpecklepyTypeError('Detector', 'optics_transmission', type(optics_transmission), 'float')

		if isinstance(quantum_efficiency, Quantity):
			self.quantum_efficiency = quantum_efficiency
		elif isinstance(quantum_efficiency, (int, float)):
			logger.warning(f"Interpreting scalar type quantum_efficiency as {quantum_efficiency} electron/ photon")
			self.quantum_efficiency = quantum_efficiency * Unit('electron / ph')
		elif isinstance(quantum_efficiency, str):
			self.quantum_efficiency = Quantity(quantum_efficiency)
		else:
			raise SpecklepyTypeError('Detector', 'quantum_efficiency', type(quantum_efficiency), 'Quantity')

		if isinstance(system_gain, Quantity):
			self.system_gain = system_gain
		elif isinstance(system_gain, (int, float)):
			logger.warning(f"Interpreting scalar type system_gain as {system_gain} electron/ ADU")
			self.system_gain = system_gain * Unit('electron / adu')
		elif isinstance(system_gain, str):
			self.system_gain = Quantity(system_gain)
		else:
			raise SpecklepyTypeError('Detector', 'system_gain', type(system_gain), 'Quantity')

		if dark_current is None or isinstance(dark_current, Quantity):
			self.dark_current = dark_current
		elif isinstance(dark_current, (int, float)):
			logger.warning(f"Interpreting scalar type dark_current as {dark_current} electron/ s")
			self.dark_current = dark_current * Unit('electron / s')
		elif isinstance(dark_current, str):
			self.dark_current = Quantity(dark_current)
		else:
			raise SpecklepyTypeError('Detector', 'dark_current', type(dark_current), 'Quantity')

		if readout_noise is None or isinstance(readout_noise, Quantity):
			self.readout_noise = readout_noise
		elif isinstance(readout_noise, (int, float)):
			logger.warning(f"Interpreting scalar type readout_noise as {readout_noise} electron")
			self.readout_noise = readout_noise * Unit('electron')
		elif isinstance(readout_noise, str):
			self.readout_noise = Quantity(readout_noise)
		else:
			raise SpecklepyTypeError('Detector', 'readout_noise', type(readout_noise), 'Quantity')

		if isinstance(saturation_level, Quantity) or saturation_level is None:
			self.saturation_level = saturation_level
		elif isinstance(saturation_level, (int, float)):
			logger.warning(f"Interpreting scalar type saturation_level as {saturation_level} electron")
			self.saturation_level = saturation_level * Unit('electron')
		elif isinstance(saturation_level, str):
			self.saturation_level = Quantity(saturation_level)
		else:
			raise SpecklepyTypeError('Detector', 'saturation_level', type(saturation_level), 'Quantity')

		# Derive secondary parameters
		self.array = np.zeros(self.shape)
		self.field_of_view = (self.shape[0] * self.pixel_scale, self.shape[1] * self.pixel_scale)

	@staticmethod
	def from_file(par_file):
		params = config.read(par_file)

		# Try known first-level keys
		for key in ['DETECTOR', 'Detector', 'detector']:
			if key in params.keys():
				return Detector(**params[key])

		# Try full parameter set
		try:
			return Detector(**params)
		except TypeError:
			raise RuntimeError(f"Could not identify parameters for initializing a Detector instance in parameter "
							   f"file {par_file}")

	@property
	def resolution(self):
		return self.pixel_scale

	@resolution.setter
	def resolution(self, value):
		self.pixel_scale = value

	def __call__(self, *args, **kwargs):
		return self.get_counts(*args, **kwargs)

	def __str__(self):
		tmp = "Detector:"
		for key in self.__dict__:
			if key == 'array':
				continue
			tmp += f"\n{key}: {self.__dict__[key]}"
		return tmp

	def get_counts(self, photon_rate, integration_time, photon_rate_resolution, debug=False):
		"""Computes the counts array from the photon rate.

		Args:
			photon_rate (Quantity):
				Passed to expose() method.
			integration_time (Quantity):
				Passed to expose() and readout() methods.
			photon_rate_resolution (Quantity):
				Angular resolution of the photon_rate array, used for resampling this to the detectors grid.
			debug (bool, optional):
				Set True for debugging. Default is False.

		Returns:
			counts (Quantity):
				Array of the shape of the detector that contains the counts measured within every pixel.
		"""
		self.expose(photon_rate=photon_rate, integration_time=integration_time,
					photon_rate_resolution=photon_rate_resolution, debug=debug)
		return self.readout(integration_time=integration_time)

	def resample(self, photon_rate, photon_rate_resolution):
		"""Resamples the photon_rate array to the angular resolution of the detector.

		Args:
			photon_rate (Quantity):
			photon_rate_resolution (Quantity):

		Returns:
			photon_rate_resampled_subset (Quantity):
				Resampled subset of the photon_rate array.
		"""

		# Assert that the photon_rate covers a larger field of view than the detector field of view
		photon_rate_field_of_view = (photon_rate.shape[0] * photon_rate_resolution,
								   photon_rate.shape[1] * photon_rate_resolution)
		if photon_rate_field_of_view[0] < self.field_of_view[0] or \
				photon_rate_field_of_view[1] < self.field_of_view[1]:
			raise ValueError(f"The field of view of the photon rate image ({photon_rate_field_of_view}) is smaller "
							 f"than that of the detector ({self.field_of_view}) in at least one dimension!")

		# Resample the photon_rate array to the detector resolution
		zoom_ratio = float(photon_rate_resolution / self.resolution)
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			photon_rate_resampled = zoom(photon_rate, zoom_ratio, order=1) * photon_rate.unit
			photon_rate_resampled = photon_rate_resampled / zoom_ratio**2  # This is necessary for flux conservation

		# Extract the central region of shape=Detector.shape
		center = (int(photon_rate_resampled.shape[0] / 2), int(photon_rate_resampled.shape[1] / 2))
		dx = int(self.shape[0] / 2)
		dy = int(self.shape[1] / 2)
		photon_rate_resampled_subset = photon_rate_resampled[center[0] - dx: center[0] + dx, center[1] - dy: center[1] + dy]

		return photon_rate_resampled_subset

	def expose(self, photon_rate, integration_time, photon_rate_resolution, debug=False):
		"""Compute the number of electrons in every pixel after the exposure.

		Args:
			photon_rate (Quantity):
				Passed to expose() method.
			integration_time (Quantity):
				Passed to expose() and readout() methods.
			photon_rate_resolution (Quantity):
				Angular resolution of the photon_rate array, used for resampling this to the detectors grid.
			debug (bool, optional):
				Set True for debugging. Default is False.

		Returns:
			electrons (Quantity):

		"""

		# Input parameters
		if isinstance(photon_rate, (int, float)):
			logger.warning(f"Interpreting scalar type photon_rate as {photon_rate} photon/ s")
			photon_rate = photon_rate * Unit('ph / s')
		elif not isinstance(photon_rate, Quantity):
			raise SpecklepyTypeError('expose', 'photon_rate', type(photon_rate), 'Quantity')

		if isinstance(integration_time, (int, float)):
			logger.warning(f"Interpreting scalar type integration_time as {integration_time} s")
			integration_time = integration_time * Unit('s')
		elif not isinstance(integration_time, Quantity):
			raise SpecklepyTypeError('expose', 'integration_time', type(integration_time), 'Quantity')

		if isinstance(photon_rate_resolution, (int, float)):
			logger.warning(f"Interpreting scalar type photon_rate_resolution as {photon_rate_resolution} arcsec")
			photon_rate_resolution = photon_rate_resolution * Unit('arcsec')
		elif not isinstance(photon_rate_resolution, Quantity):
			raise SpecklepyTypeError('expose', 'photon_rate_resolution', type(photon_rate_resolution), 'Quantity')

		# Resample the photon rate to the detector resolution
		photon_rate = self.resample(photon_rate=photon_rate, photon_rate_resolution=photon_rate_resolution)
		photons = photon_rate * integration_time
		if debug:
			imshow(photons, title='photons')

		# Compute photon shot noise with Poisson statistics
		photons = np.random.poisson(photons.value) * photons.unit

		# Incorporate efficiencies
		if self.optics_transmission is not None:
			photons = photons * self.optics_transmission
		electrons = photons * self.quantum_efficiency
		if debug:
			imshow(electrons, title='electrons')

		# Limit to the saturation level of the detector
		if self.saturation_level is not None:
			electrons = np.minimum(electrons, self.saturation_level) # * self.system_gain)
		electrons = np.round(electrons)
		self.array = electrons
		return electrons

	def readout(self, integration_time, reset=True):
		"""Computes the readout of the detector and returns the ADUs for every pixel.

		Args:
			integration_time (Quantity):
				Integration time of an individual exposure.
			reset (bool, optional):
				If set to True, then the electron count of every pixel is reset to zero for the next exposure. Default
				is True.

		Returns:
			counts (Quantity):
				Array of the ADUs for every pixel.
		"""

		# Read copy and clear the array
		electrons = self.array

		# Apply dark_current and readout noise following Poisson or Gaussian statistics
		if self.dark_current is not None:
			dark_current_electrons = np.round(np.random.poisson(self.dark_current.value, self.shape)) \
									 * self.dark_current.unit * integration_time
			electrons = electrons + dark_current_electrons
		if self.readout_noise is not None:
			readout_electrons = np.round(np.random.normal(0.0, self.readout_noise.value, self.shape)) \
								* self.readout_noise.unit
			electrons = electrons + readout_electrons

		# Convert into ADU
		counts = electrons / self.system_gain
		if self.saturation_level is not None:
			return np.minimum(counts, self.saturation_level / self.system_gain)

		# Reset the detector array
		if reset:
			self.array = np.zeros(self.shape)

		return counts.decompose()
