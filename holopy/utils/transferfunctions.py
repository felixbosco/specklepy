import numpy as np
from numpy.fft import fft2, fftshift

def otf(psf, inverse=False):
    if inverse:
        return fftshift(ifft2(psf))
    else:
        return fftshift(fft2(psf))

def mtf(psf):
    return np.abs(otf(psf))

def powerspec(array):
    return np.square(mtf(array))

def psf(aperture):
    return powerspec(aperture)

def ft(arg):
    return otf(arg)

def powerspec1d(array, flatten=True, average=True, pixel_scale=None):
    center = (array.shape[0] / 2, array.shape[1] / 2)
    xx, yy = np.mgrid[:array.shape[0], :array.shape[1]]
    Fourier_radius = np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))
    Fourier_radius = np.ma.masked_greater(Fourier_radius, center[0])
    signal = powerspec(array)
    signal = np.ma.masked_array(signal, mask=Fourier_radius.mask)

    # Flatten the arrays
    if flatten:
        Fourier_radius = Fourier_radius.reshape((-1))
        signal = signal.reshape((-1))

        # Average signal along Fourier distance
        if average:
            xdata = np.unique(Fourier_radius)
            ydata = np.zeros(xdata.shape)
            for index, value in enumerate(xdata):
                ydata[index] = np.mean(signal[np.where(Fourier_radius == value)])
        else:
            xdata = Fourier_radius
            ydata = signal

    # Apply pixel_scale
    if pixel_scale is not None:
        size = array.shape[0]
        freq1 = np.fft.fftfreq(size, pixel_scale.value)[1]
        xdata = xdata * freq1 / pixel_scale.unit
    else:
        pass

    return xdata, ydata
