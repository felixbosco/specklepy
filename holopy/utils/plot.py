import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from holopy.config import plotting
from holopy.utils import transferfunctions as tf


def imshow(image, title=None, norm=None):
    if norm == 'log':
        norm = LogNorm()
    plt.figure()
    plt.imshow(image, norm=norm)
    plt.title(title)
    plt.colorbar(pad=0.0)
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


def plot_simple(xdata, ydata, title=None, xlabel=None, ylabel=None):
    plt.plot(xdata, ydata)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    plt.close()
