import numpy as np
import matplotlib.pyplot as plt
from holopy.config import plotting
from holopy.utils import transferfunctions as tf

def imshow(image, title=None, norm=None):
    plt.figure()
    plt.imshow(image, norm=norm)
    plt.title(title)
    plt.colorbar(pad=0.0)
    plt.show()
    plt.close()

def plot_powerspec1d(image, average=True):
    plt.figure()
    if len(image.shape) is 2:
        xdata, ydata = tf.powerspec1d(image, average=average)
        plt.plot(xdata, ydata, '-')
    elif len(image.shape) is 3:
        for frame in image:
            xdata, ydata = tf.powerspec1d(frame, average=average)
            plt.plot(xdata, ydata, '-')
    plt.xlabel('uv Radius (pixel)')
    plt.ylabel('Signal')
    # plt.ylim(0.0)
    plt.grid()
    plt.show()
    plt.close()
