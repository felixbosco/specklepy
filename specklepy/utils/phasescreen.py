import numpy as np
from scipy.fftpack import fft2, fftshift


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

    def generate_screen(self, size=None):
        if size is not None:
            self.size = size
        self.complex_screen = fft2(fftshift(self.psd()))
        return self.complex_screen

    def generate(self, number=None, size=None):
        if number is None or number == 1:
            return self.generate(size=size).real
        else:
            screens = []
            for n in range(number):
                screens.append(self.generate_screen(size=size).real)
            return screens


class Atmosphere(object):

    def __init__(self, screens, speeds, fractions=None):

        if len(screens) != len(speeds):
            raise ValueError(f"The number of phase screens ({len(screens)}) has to be the same as number of wind "
                             f"speeds ({len(speeds)})!")
        self.screens = np.array(screens)
        self.speeds = speeds
        self.n_layers = len(screens)

        if fractions is None:
            self.screen_weights = fractions
        else:
            self.screen_weights = np.expand_dims(np.array(fractions), axis=(1, 2))

        # Initialize secondary parameters
        self.step = 0

    def __repr__(self):
        return f"Atmosphere(screens={self.n_layers}, step={self.step})"

    @property
    def screen_shape(self):
        return self.screens[0].shape

    @property
    def weighted_screens(self):
        if self.screen_weights is None:
            return self.screens
        else:
            return np.multiply(self.screen_weights, self.screens)

    def integrate_layers(self):
        return np.sum(self.weighted_screens, axis=0)

    def evolve(self):
        for n in range(self.n_layers):
            self.screens[n] = np.roll(self.screens[n], shift=self.speeds[n])
        self.step += 1


class PSFIterator(object):

    def __init__(self, radius, atmosphere):
        # Store input
        self.radius = radius
        self.atmosphere = atmosphere

        # Initialize secondary parameters
        self.aperture = self.initialize_circular_aperture()

    def __repr__(self):
        return f"PSFIterator()"

    @property
    def size(self):
        return self.atmosphere.screen_shape[0]

    def initialize_circular_aperture(self):
        xx, yy = np.mgrid[:self.size, :self.size]
        return np.square(xx - self.size / 2) + np.square(yy - self.size / 2) <= np.square(self.radius)

    @property
    def complex_aperture(self):
        return np.multiply(self.aperture, np.exp(1j * self.atmosphere.integrate_layers()))

    @property
    def psf(self):
        return fftshift(np.square(np.abs(fft2(fftshift(self.complex_aperture)))))

    def next(self):
        self.atmosphere.evolve()
        return self.psf
