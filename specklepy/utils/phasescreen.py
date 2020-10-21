import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift

from specklepy.utils.transferfunctions import psf


class PhaseScreen(object):

    def __init__(self, k0=3, norm=3, size=256):
        self.k0 = k0
        self.norm = norm
        self.size = size
        self.complex_screen = None

    @property
    def shape(self):
        return self.size, self.size

    @property
    def screen(self):
        return self.complex_screen.real

    def random_phase(self):
        return np.random.rand(self.size, self.size) * 2 * np.pi

    def amplitude(self):
        return np.power(np.square(self.initialize_radii() + np.square(self.k0)), -11 / 12)

    def initialize_radii(self, center=None):
        if center is None:
            center = self.size / 2
        if isinstance(center, (int, float)):
            center = tuple([center, center])
        xx, yy = np.mgrid[:self.size, :self.size]
        return np.sqrt(np.square(xx - center[1]) + np.square(yy - center[0]))

    def psd(self):
        return self.norm * self.amplitude() * np.exp(1j * self.random_phase())

    def generate(self, size=None):
        if size is not None:
            self.size = size
        self.complex_screen = fft2(fftshift(self.psd()))

    def generate_n(self, number_screens, size=None):
        screens = []
        for n in range(number_screens):
            self.generate(size=size)
            screens.append(self.screen)
        return screens


class PSFIterator(object):

    def __init__(self, width, screens, speeds, fractions=None):
        # Store input
        self.width = width

        if len(screens) != len(speeds):
            raise ValueError(f"The number of phase screens ({len(screens)}) has to be the same as number of wind "
                             f"speeds ({len(speeds)})!")
        else:
            self.screens = screens
            self.speeds = speeds
            self.n_layers = len(screens)

        if fractions is None:
            self.fractions = self.initialize_fractions()
        else:
            self.fractions = fractions

        # Initialize secondary parameters
        self.indexes = self.initialize_indexes()
        self.aperture = self.initialize_circular_aperture()
        self.step = 0

    def __repr__(self):
        return f"PSFIterator(screens={self.n_layers}, step={self.step})"

    @property
    def radius(self):
        return self.width // 2

    def initialize_fractions(self):
        return [1 / self.n_layers] * self.n_layers

    def initialize_indexes(self):
        return [(self.width // 2, self.width // 2)] * self.n_layers

    def initialize_circular_aperture(self):
        xx, yy = np.mgrid[:self.width, :self.width]
        rr = np.sqrt(np.square(xx - self.radius) + np.square(yy - self.radius))
        aperture = (rr <= self.radius).astype(float)
        return aperture #/ np.sum(aperture)

    def grab_layer(self, index):
        x0, y0 = self.indexes[index]
        return self.screens[index][x0 - self.radius: x0 + self.radius, y0 - self.radius: y0 + self.radius]

    def integrate_layers(self):
        phase = self.grab_layer(0)
        for layer in range(1, self.n_layers):
            phase += self.grab_layer(layer) * self.fractions[layer]
        return phase

    def complex_aperture(self):
        # return np.multiply(self.aperture, np.exp(1j * self.integrate_layers()))
        return np.add(self.aperture, 1j * self.integrate_layers())

    def psf(self):
        return np.square(np.abs(fftshift(ifft2(self.complex_aperture()))))
