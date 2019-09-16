import numpy as np
from numpy.fft import fft2, fftshift
# from holopy.utils.plot import imshow

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

def powerspec1d(array, flatten=True, average=True, plot_intermediate=False):
    center = (array.shape[0] / 2, array.shape[1] / 2)
    xx, yy = np.mgrid[:array.shape[0], :array.shape[1]]
    Fourier_distance = np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))
    # if plot_intermediate:
    #     imshow(Fourier_distance, title='Fourier distance')
    signal = powerspec(array)
    # if plot_intermediate:
    #     imshow(signal, title='Signal')

    # Flatten the arrays
    if flatten:
        Fourier_distance = Fourier_distance.reshape((-1))
        signal = signal.reshape((-1))

        # Average signal along Fourier distance
        if average:
            xdata = np.unique(Fourier_distance)
            ydata = np.zeros(xdata.shape)
            for index, value in enumerate(xdata):
                ydata[index] = np.mean(signal[np.where(Fourier_distance == value)])
            return xdata, ydata
    return Fourier_distance, signal
