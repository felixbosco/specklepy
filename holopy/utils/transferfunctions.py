import numpy as np
from numpy.fft import fft2, fftshift
from holopy.utils.plot import imshow

def otf(psf):
    return fftshift(fft2(psf))

def mtf(psf):
    return np.abs(otf(psf))

def powerspec(array):
    return np.square(mtf(array))

def psf(aperture):
    return powerspec(aperture)

def ft(arg):
    return otf(arg)

def powerspec1d(array, plot_intermediate=False):
    center = (array.shape[0] / 2, array.shape[1] / 2)
    print(array.shape)
    print(center)
    xx, yy = np.mgrid[:array.shape[0], :array.shape[1]]
    Fourier_distance = np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))
    if plot_intermediate:
        imshow(Fourier_distance, title='Fourier distance')
    signal = powerspec(array)
    if plot_intermediate:
        imshow(signal, title='Signal')
    return Fourier_distance.reshape((-1)), signal.reshape((-1))
