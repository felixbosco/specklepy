from IPython import embed
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import sys

from astropy.io import fits
from astropy.io.registry import IORegistryError
from astropy.table import Table

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.plotting.plots import save_figure


class Plot(object):

    def __init__(self, x_data=None, y_data=None, image_data=None, file=None, layout=None, debug=False, **kwargs):

        # Store data
        self.x_data = x_data
        self.y_data = y_data
        self.image_data = image_data
        self.file = file
        self.layout = layout
        self.debug = debug

        # Create default attributes
        self.colorbar = None

        # Create figure
        self.figure, self.axes = plt.subplots(**kwargs)
        if not isinstance(self.axes, list):
            self.axes = [self.axes]
        self.plot_data()
        self.plot_image()

    @classmethod
    def from_file(cls, file_name, extension=None, columns=None, format=None, layout=None, debug=False):
        # Does file exist?
        if not os.path.exists(path=file_name):
            sys.tracebacklimit = 0
            raise FileNotFoundError(f"File {file_name} not found!")

        # Set data defaults
        x_data = None
        y_data = None
        image_data = None
        file = cls.default_file_name(input=file_name)

        # Apply fall-back value for format
        if format is None:
            format = 'ascii.fixed_width'

        # Identify file type and read data
        root, ext = os.path.splitext(file_name)
        try:
            if 'fits' in ext:
                table = Table.read(file_name, hdu=extension)
            else:
                table = Table.read(file_name, format=format)
            if columns is None:
                columns = table.colnames[0:2]
            x_data = table[columns[0]].data
            y_data = table[columns[1]].data
        except IORegistryError as e:
            sys.tracebacklimit = 0
            raise e
        except ValueError:
            image_data = fits.getdata(file_name, extension)

        return cls(x_data=x_data, y_data=y_data, image_data=image_data, file=file, layout=layout, debug=debug)

    def set_height(self, val):
        if isinstance(val, str):
            if val == 'text':
                val = 10
            elif val == 'column':
                val = 5
            else:
                raise SpecklepyValueError('Plot.set_height', argname='val', argvalue=val, expected="'text' or 'column'")
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

    def plot_image(self, axis=0, **colorbar_kwargs):
        if self.image_data is not None:
            ax = self.axes[axis]
            img = ax.imshow(self.image_data, origin='lower')

            # Create color bar
            if 'pad' in colorbar_kwargs:
                pad = colorbar_kwargs['pad']
            else:
                pad = 0.0
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=pad)
            self.colorbar = plt.colorbar(img, cax=cax, ax=ax, pad=pad)
            if 'label' in colorbar_kwargs:
                self.colorbar.set_label(colorbar_kwargs['label'])

    def apply_layout(self, layout=None):
        if layout is None:
            return
        else:
            self.set_width(val=layout)
            self.set_height(val=layout)
            plt.rc('font', family='serif', size=10)

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
