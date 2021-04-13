import numpy as np


class LinearGrid(object):

    def __init__(self, min, max, bin_size):
        self.bins = np.arange(start=min, stop=max+bin_size, step=bin_size)

    def __len__(self):
        return len(self.bins) - 1

    @classmethod
    def from_array(cls, array, bin_size):
        finite = np.isfinite(array)
        a_max = np.max(array[finite])
        a_min = np.min(array[finite])
        return cls(min=a_min, max=a_max, bin_size=bin_size)

    def average(self, x_data, y_data, e_data=None):
        x_out = np.empty((len(self)))
        y_out = np.empty((len(self)))
        e_out = np.empty((len(self)))

        for b, bin_lower in enumerate(self.bins[:-1]):
            bin_upper = self.bins[b + 1]
            indexes = np.where(np.logical_and(bin_lower <= x_data, x_data < bin_upper))

            x_out[b] = np.mean(x_data[indexes])
            y_out[b] = np.mean(y_data[indexes])
            if e_data is None:
                e_out[b] = np.std(y_data[indexes])
            else:
                e_out[b] = np.sqrt(np.sum(np.square(e_data[indexes])))

        return x_out, y_out, e_out

    @classmethod
    def regrid(cls, bin_size, x_data, y_data, e_data=None):
        obj = cls.from_array(x_data, bin_size=bin_size)
        return obj.average(x_data=x_data, y_data=y_data, e_data=e_data)
