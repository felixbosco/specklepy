import numpy as np

import matplotlib.pyplot as plt

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError


class Plot(object):

    def __init__(self, *data, file=None, **kwargs):

        self.figure = plt.figure(**kwargs)
        self.data = data
        self.file = file

        # Digest input data
        self.data_type = None
        if isinstance(data, list):
            pass
        elif isinstance(data, np.ndarray):
            if data.ndim == 2:
                self.data_type = 'image'
        else:
            raise SpecklepyTypeError('Plot', argname='data', argtype=type(data), expected="'list' or 'np.array'")

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

    def save(self, file=None):
        if file is None:
            file = self.file
        self.figure.savefig(file)

    def show(self):
        self.figure.show()