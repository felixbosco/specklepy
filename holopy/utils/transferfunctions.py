from numpy import abs, square
from numpy.fft import fft2, fftshift

def otf(psf):
    return fftshift(fft2(psf))

def mtf(psf):
    return abs(otf(psf))

def powerspec(array):
    return square(mtf(array))

def psf(aperture):
    return powerspec(aperture)

def ft(arg):
    return otf(arg)

def powerspec1d(array):
    center = (array.shape[0] / 2, array.shape[1])
    xx, yy = np.mgrid[:array.shape[0], :array.shape[1]]
    Fourier_distance = np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))
    signal = powerspec(array)
    return Fourier_distance.reshape((-1)), signal.reshape((-1))
