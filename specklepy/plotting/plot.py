from IPython import embed
import numpy as np
import os
import sys

import matplotlib.pyplot as plt

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.table import read_table
from specklepy.plotting.plots import save_figure


class Plot(object):

    def __init__(self, x_data=None, y_data=None, image_data=None, file=None, debug=False, **kwargs):

        # Store data
        self.x_data = x_data
        self.y_data = y_data
        self.image_data = image_data
        self.file = file
        self.debug = debug

        # Create figure
        self.figure, self.axes = plt.subplots(**kwargs)
        if not isinstance(self.axes, list):
            self.axes = [self.axes]
        self.plot_data()
        self.plot_image()

    @classmethod
    def from_file(cls, file_name, extension=None, columns=None, debug=False):
        # Does file exist?
        if not os.path.exists(path=file_name):
            sys.tracebacklimit = 0
            raise FileNotFoundError(f"File {file_name} not found!")

        # Set data defaults
        x_data = None
        y_data = None
        image_data = None
        file = cls.default_file_name(input=file_name)

        # Identify file type and read data
        root, ext = os.path.splitext(file_name)
        if 'fits' in ext:
            pass
        else:
            table = read_table(file_name=file_name)
            if columns is None:
                columns = table.colnames[0:2]
            x_data = table[columns[0]].data
            y_data = table[columns[1]].data

        return cls(x_data=x_data, y_data=y_data, image_data=image_data, file=file, debug=debug)

    def set_height(self, val):
        self.figure.set_figheight(val=val)

    def set_width(self, val):
        if isinstance(val, str):
            if val == 'text':
                val = 10
            elif val == 'column':
                val = 5
            else:
                raise SpecklepyValueError('Plot.set_width', argname='val', argvalue=val, expected="'text' or 'column'")
        elif not isinstance(val, (int, float)):
            raise SpecklepyTypeError('Plot.set_width', argname='val', argtype=type(val), expected='str or float')
        self.figure.set_figwidth(val=val)

    def plot_data(self, axis=0):
        if self.y_data is not None:
            self.axes[axis].plot(self.x_data, self.y_data)

    def plot_image(self, axis=0):
        if self.image_data is not None:
            self.axes[axis].imshow(self.image_data)

    def save(self, file=None):
        if file is None:
            file = self.file
        save_figure(file)

    @staticmethod
    def show():
        plt.show()

    @staticmethod
    def default_file_name(input=None):
        if input is None:
            return 'plot.png'
        else:
            root, ext = os.path.splitext(os.path.basename(input))
            return "plot_" + root + '.png'
