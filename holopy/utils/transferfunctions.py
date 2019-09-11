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
