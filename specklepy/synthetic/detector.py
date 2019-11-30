import sys
import numpy as np
from astropy import units as u
from scipy.ndimage import zoom
import warnings

class Detector(object):

	"""Class carrying information of an astronomical detector.

	Attributes:
		shape (tuple, dtype=int): pixel times pixel
		pixel_size (astropy.units.Quantity): in arcsec or arcsec per pix
		quantum_efficiency ():
		system_gain (astropy.units.Quantity):
		
	Optional attributes are:
		readout_noise (astropy.units.Quantity):
		dark_current (astropy.units.Quantity):
		saturation_level (astropy.units.Quantity):
		optics_transmission (astropy.units.Quantity):

	Future features:
		Method 'window': Shall enable the interactive windowing of a detector.
		Attribute dictionary 'readout_modes': Shall enable the flexible use of different
			readout modes with different parameters.
	"""

	__name__ = 'detector'

	def __init__(self, shape, pixel_size, quantum_efficiency=1.*u.electron/u.ph, system_gain=1.*u.electron/u.adu, **kwargs):
		# Read input parameters
		self.shape = shape
		self.pixel_size = pixel_size
		self.quantum_efficiency = quantum_efficiency
		self.system_gain = system_gain
		for key in kwargs:
			self.__setattr__(key, kwargs[key])

		# Compute secondary parameters
		self.array = np.zeros(self.shape)
		self.FoV = (self.shape[0] * self.pixel_size, self.shape[1] * self.pixel_size)


	@property
	def resolution(self):
		return self.pixel_size


	@resolution.setter
	def resolution(self, value):
		self.pixel_size = value


	def __call__(self, photon_rate_density_array, integration_time, target_FoV, compute_photon_shot_noise=True):
		self.expose(photon_rate_density_array, integration_time, target_FoV, compute_photon_shot_noise=compute_photon_shot_noise)
		return self.readout(integration_time=integration_time)


	def __str__(self):
		tmp = "Detector:\n"
		for key in self.__dict__:
			if key == 'array':
				continue
			tmp += "{}: {}\n".format(key, self.__dict__[key])
		return tmp


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


	def expose(self, photon_rate_density_array, integration_time, target_FoV, compute_photon_shot_noise=True, verbose=0):
		tmp = self.resample(photon_rate_density_array, target_FoV) * integration_time
		tmp *= self.quantum_efficiency
		if hasattr(self, 'optics_transmission'):
			tmp *= self.optics_transmission
		if compute_photon_shot_noise:
			try:
				tmp = np.random.poisson(tmp.value) * tmp.unit
			except ValueError as e:
				if verbose > 0:
					print('Bypassed ValueError ({}) in np.random.poisson() by substituting values smaller than zero by zero.'.format(e))
				tmp = np.random.poisson(np.maximum(tmp.value, 0.0)) * tmp.unit
		if hasattr(self, 'saturation_level'):
			tmp = np.minimum(tmp, self.saturation_level * self.system_gain)
		self.array = np.round(tmp)


	def readout(self, integration_time, reset=True):
        # Read copy and clear the array
		tmp = self.array
		if hasattr(self, 'dark_current'):
			tmp += np.random.poisson(self.dark_current.value, self.shape) * self.dark_current.unit * u.pix * integration_time
		if hasattr(self, 'readout_noise'):
			tmp += np.round(np.random.normal(0.0, self.readout_noise.value, self.shape) ) * self.readout_noise.unit * u.pix
		tmp /= self.system_gain
		if hasattr(self, 'saturation_level'):
			return np.minimum(tmp, self.saturation_level)
		# Reset the detector
		if reset:
			self.array = np.zeros(self.shape)
		return tmp.decompose()


	def window(self):
		pass
