import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs

from specklepy.logging import logging
from specklepy.utils import transferfunctions as tf



# Define plotting defaults
font = {
    'family' : 'serif',
    #'weight' : 'bold',
    #'size'   : 'larger'
    }

plt.rc('font', **font)  # pass in the font dict as kwargs



def imshow(image, title=None, norm=None):
    if norm == 'log':
        norm = clrs.LogNorm()
    plt.figure()
    plt.imshow(image, norm=norm)
    plt.title(title)
    plt.colorbar(pad=0.0)
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
        logging.info("Interpreting color tuple () as RGB values.")
        hsv_color = clrs.rgb_to_hsv(color)
    else:
        raise TypeError("The function desaturate_color received color argument \
                        of type {}, but needs to be str!".format(type(color)))

    if not isinstance(ncolors, int):
        raise TypeError("The function desaturate_color received ncolor argument\
                        of type {}, but needs to be int!".format(type(ncolors)))

    if not isinstance(saturation_min, float):
        raise TypeError("The function desaturate_color received saturation_min\
                        argument of type {}, but needs to be float!".format(type(saturation_min)))

    if saturation_values is None:
        saturation_values = np.linspace(hsv_color[1], saturation_min, num=ncolors)
    elif isinstance(saturation_values, list):
        pass # list is correct, nothing to adapt
    elif isinstance(saturation_values, float):
        saturation_values = list(saturation_values)
    else:
        raise TypeError("The function desaturate_color received \
                        saturation_values argument of type {}, but needs to be list!".format(type(saturation_values)))


    # Create list of colors with varied saturation values
    colors = []
    for saturation_value in saturation_values:
        color = clrs.to_rgb((hsv_color[0], saturation_value, hsv_color[2]))
        colors.append(color)

    return colors
