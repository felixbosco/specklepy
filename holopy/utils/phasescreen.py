import numpy as np
from scipy.fftpack import fft2, fftshift

class PhaseScreen(object):

    def __init__(self, k0=3, norm=3, **kwargs):
        self.k0 = k0
        self.norm = norm

    @property
    def amplitude(self):
        return self.amp

    def __call__(self):
        return self.norm * self.amp * np.exp(1j * self.phase)

    def make_distance_map(self, size, radius=None):
        if radius is None:
            radius = size / 2
        xx, yy = np.mgrid[:size, :size]
        return np.sqrt(np.square(xx - radius) + np.square(yy - radius))

    def make_aperture(self, size, radius, amplitude=1):
        return (self.make_distance_map(size=size, radius=radius) <= radius)

    def make_phase(self, size):
        #yy, xx = np.ogrid[-size: size+1, -size: size+1]
        #kr = np.sqrt(xx**2 + yy**2)
        kr = self.make_distance_map(size=size)
        psd_amp = np.power(np.square(kr) + np.square(self.k0), -11 / 12)
        psd_pha = np.random.rand(size, size) * 2 * np.pi
        psd = self.norm * psd_amp * np.exp(1j * psd_pha)
        self.screen = fft2(fftshift(psd)).real
        return self.screen
