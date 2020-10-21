import numpy as np
from scipy.fftpack import fft2, fftshift


class PhaseScreen(object):

    def __init__(self, k0=3, norm=3, **kwargs):
        self.k0 = k0
        self.norm = norm

    @staticmethod
    def random_phase(size=256):
        return np.random.rand(size, size) * 2 * np.pi

    def amplitude(self, size=256):
        return np.power(np.square(self.initialize_radii(size=size) + np.square(self.k0)), -11 / 12)

    @staticmethod
    def initialize_radii(size=256):
        center = size / 2
        xx, yy = np.mgrid[:size, :size]
        return np.sqrt(np.square(xx - center) + np.square(yy - center))

    def psd(self, size=256):
        return self.norm * self.amplitude(size=size) * np.exp(1j * self.random_phase(size=size))

    def generate(self, size=256):
        return fft2(fftshift(self.psd(size=size))).real
