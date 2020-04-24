import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import astropy.units as u

from specklepy.logging import logger
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.utils import transferfunctions as tf


# Define plotting defaults
plt.rc('figure', figsize=(11.69, 8.27))  # A4 format in inches
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
plt.rc('font', family='serif', size=14)
DPI = 400


def imshow(image, title=None, norm=None, colorbar_label=None, saveto=None, maximize=False):
    """Shows a 2D image.

    Args:
        image (np.ndarray, ndim=2):
            Image to be plotted.
        title (str, optional):
            Plot title. Default is None.
        norm (str, optional):
            Can be set to 'log', for plotting in logarithmic scale. Default is
            None.
        colorbar_label (str, optional):
            Label of the color bar. Default is None.
        saveto (str, optional):
            Path to save the plot to. Default is None.
        maximize (bool, optional):
            Set true for showing the plot on full screen. Default is False.
    """

    if isinstance(image, np.ndarray):
        if image.ndim != 2:
            raise SpecklepyValueError('imshow()', 'image.ndim', image.ndim, '2')
        if isinstance(image, u.Quantity):
            unit = image.unit
            colorbar_label = "({})".format(unit)
            image = image.value
    else:
        raise SpecklepyTypeError('imshow()', 'image', type(image), 'np.ndarray')

    if norm == 'log':
        norm = clrs.LogNorm()
    plt.figure()
    plt.imshow(image, norm=norm)
    plt.title(title)
    if maximize:
        maximize_plot()

    # Colorbar
    cbar = plt.colorbar(pad=0.0)
    if colorbar_label is not None:
        cbar.set_label(colorbar_label)

    if saveto is not None:
        plt.savefig(saveto, dpi=300)

    plt.show()
    plt.close()


def plot_simple(xdata, ydata, title=None, xlabel=None, ylabel=None):
    plt.plot(xdata, ydata)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    plt.close()


def plot_powerspec1d(image, title=None, average=True, pixel_scale=None):
    plt.figure()
    if image.ndim is 2:
        xdata, ydata = tf.powerspec1d(image, average=average, pixel_scale=pixel_scale)
        plt.plot(xdata, ydata, '-')
    elif image.ndim is 3:
        for frame in image:
            xdata, ydata = tf.powerspec1d(frame, average=average, pixel_scale=pixel_scale)
            plt.plot(xdata, ydata, '-')
    plt.xlabel('uv Radius (pixel)')
    plt.ylabel('Signal')
    plt.title(title)
    plt.grid()
    plt.show()
    plt.close()


def maximize_plot():
    mng = plt.get_current_fig_manager()
    backend = mpl.get_backend()
    if backend == 'Qt4Agg' or backend == 'Qt5Agg':
        mng.window.showMaximized()
    elif backend == 'wxAgg':
        mng.frame.Maximize(True)
    elif backend == 'TkAgg':
        try:
            mng.window.state('zoomed')
        except:
            logger.warning("Could not maximize plot")
    else:
        raise RuntimeWarning("Maximizing plot is not possible with matplotlib backend {}".format(backend))


def desaturate_color(color, ncolors=1, saturation_values=None, saturation_min=0.1):
    """Desaturates a color and returns a list of desaturated colors.

    Args:
        color (str or ...):
        ncolors (int): Number of returned colors.
        saturation_values (None or list, dtype=float):
        saturation_min (float): Minimum value of saturation.

    Returns:
        colors (list): List of RGB represnations of colors with length ncolors.
    """

    # Input parameters
    if isinstance(color, str):
        rgb_color = clrs.to_rgb(color)
        hsv_color = clrs.rgb_to_hsv(rgb_color)
    elif isinstance(color, tuple):
        logger.info("Interpreting color tuple () as RGB values.")
        hsv_color = clrs.rgb_to_hsv(color)
    else:
        raise SpecklepyTypeError('desaturate_color()', 'color', type(color), 'str')

    if not isinstance(ncolors, int):
        raise SpecklepyTypeError('desaturate_color()', 'ncolors', type(ncolors), 'int')

    if not isinstance(saturation_min, float):
        raise SpecklepyTypeError('desaturate_color()', 'saturation_min', type(saturation_min), 'float')

    if saturation_values is None:
        saturation_values = np.linspace(hsv_color[1], saturation_min, num=ncolors)
    elif isinstance(saturation_values, list):
        pass  # list is correct, nothing to adapt
    elif isinstance(saturation_values, float):
        saturation_values = list(saturation_values)
    else:
        raise SpecklepyTypeError('desaturate_color()', 'saturation_values', type(saturation_values), 'list')

    # Create list of colors with varied saturation values
    colors = []
    for saturation_value in saturation_values:
        color = clrs.hsv_to_rgb((hsv_color[0], saturation_value, hsv_color[2]))
        colors.append(color)

    return colors


def psf_profile_plot(files, normalize=None, maximize=False):
    """Plots the psf profile data from a file.

    Args:
        files (str):
            Name of the file to extract the data from.
        normalize (str, optional):
            Normalization mode, can be None, 'peak' and 'last'. Default is None.
        maximize (bool, optional:
            Show plots on full screen. Default is False.
    """

    # Input parameters
    if isinstance(files, str):
        plt.title(files)
        files = [files]
    if not isinstance(files, list):
        raise SpecklepyTypeError('psf_profile_plot()', 'files', type(files), 'list')

    for file in files:
        xdata, ydata = np.loadtxt(file).transpose()

        plt.xlabel("Radius (pix)")
        if normalize == 'peak':
            plt.ylabel("Peak normalized flux")
            ydata /= ydata[0]
        elif normalize == 'last':
            plt.ylabel("Sky normalized flux")
            ydata /= ydata[-1]
        else:
            plt.ylabel("Flux")

        plt.plot(xdata, ydata)
    if maximize:
        maximize_plot()
    plt.show()
    plt.close()


def encircled_energy_plot(files, normalize=None, maximize=False):
    """Plots the encircled energy data from a file.

    Args:
        files (str):
            Name of the file to extract the data from.
        normalize (str, optional):
            Normalization mode, can be None, 'peak' and 'last'. Default is None.
        maximize (bool, optional:
            Show plots on full screen. Default is False.
    """

    # Input parameters
    if isinstance(files, str):
        plt.title(files)
        files = [files]
    if not isinstance(files, list):
        raise SpecklepyTypeError('encircled_energy_plot()', 'files', type(files), 'list')

    for file in files:
        xdata, ydata = np.loadtxt(file).transpose()

        plt.xlabel("Radius (pix)")
        if normalize == 'peak':
            plt.ylabel("Peak normalized flux")
            ydata /= ydata[0]
        elif normalize == 'last':
            plt.ylabel("Sky normalized flux")
            ydata /= ydata[-1]
        else:
            plt.ylabel("Flux")

        plt.plot(xdata, ydata)
    if maximize:
        maximize_plot()
    plt.show()
    plt.close()
