import numpy as np

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.logging import logger


class Box(object):

    def __init__(self, indexes=None):

        # Initialize attributes with kwargs
        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None

        # Interpret indexes input and overwrite
        if indexes is not None:
            if not isinstance(indexes, (list, np.ndarray)):
                raise SpecklepyTypeError('Box', argname='*args', argtype=type(indexes), expected='list')
            else:
                if len(indexes) == 1:
                    if len(indexes[0]) == 1:
                        raise SpecklepyValueError
                    else:
                        logger.debug("Using the same limits for both axes")
                        self.x_min = indexes[0][0]
                        self.x_max = indexes[0][1]
                        self.y_min = indexes[0][0]
                        self.y_max = indexes[0][1]
                elif len(indexes) == 2:
                    if isinstance(indexes[0], int) and isinstance(indexes[1], int):
                        logger.debug("Using the same limits for both axes")
                        self.x_min = indexes[0]
                        self.x_max = indexes[1]
                        self.y_min = indexes[0]
                        self.y_max = indexes[1]
                    elif isinstance(indexes[0], list) and isinstance(indexes[1], list):
                        try:
                            self.x_min = indexes[0][0]
                            self.x_max = indexes[0][1]
                        except IndexError:
                            pass
                        try:
                            self.y_min = indexes[1][0]
                            self.y_max = indexes[1][1]
                        except IndexError:
                            pass
                elif len(indexes) == 4:
                    self.x_min = indexes[0]
                    self.x_max = indexes[1]
                    self.y_min = indexes[2]
                    self.y_max = indexes[3]

        # Assert integer indexes
        self.round_indexes()

    def round_indexes(self):
        for attr in ['x_min', 'x_max', 'y_min', 'y_max']:
            try:
                self.__setattr__(attr, int(round(self.__getattribute__(attr))))
            except TypeError:
                pass

    def __repr__(self):
        return f"[[{self.x_min}, {self.x_max}], [{self.y_min}, {self.y_max}]]"

    def __call__(self, array):  #, update=False):

        # Abort if Box is undefined
        if self.x_min is None and self.x_max is None and self.y_min is None and self.y_max is None:

            return array

        else:

            return array[self.y_min: self.y_max, self.x_min: self.x_max]

    def __copy__(self):
        return Box([self.x_min, self.x_max, self.y_min, self.y_max])

    @property
    def shape(self):
        if None in [self.x_max, self.x_min, self.y_max, self.y_min]:
            return None
        else:
            return self.y_max - self.y_min, self.x_max - self.x_min

    @classmethod
    def centered_at(cls, x0, y0, radius):
        return cls([x0 - radius, x0 + radius + 1, y0 - radius, y0 + radius + 1])

    def shift(self, shift, copy=True):
        # Create a copy if requested
        if copy:
            obj = self.__copy__()
        else:
            obj = self

        # Shift entries
        if shift is not None:
            obj.x_min += shift[1]
            obj.x_max += shift[1]
            obj.y_min += shift[0]
            obj.y_max += shift[0]

        obj.round_indexes()

        return obj

    def crop_to_shape(self, shape):
        if self.x_min < 0:
            self.x_min = None
        if self.y_min < 0:
            self.y_min = None
        if self.x_max >= shape[1]:
            self.x_max = None
        if self.y_max >= shape[0]:
            self.y_max = None

    def transpose(self):
        x_min = self.x_min
        x_max = self.x_max
        self.x_min = self.y_min
        self.x_max = self.y_max
        self.y_min = x_min
        self.y_max = x_max


class Box3D(object):

    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max

    def round_indexes(self):
        for attr in ['x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max']:
            try:
                self.__setattr__(attr, int(round(self.__getattribute__(attr))))
            except TypeError:
                pass

    def __repr__(self):
        return f"[[{self.x_min}, {self.x_max}], [{self.y_min}, {self.y_max}], [{self.z_min}, {self.z_max}]]"

    def __call__(self, array):  #, update=False):

        # Abort if Box is undefined
        if self.x_min is None and self.x_max is None and self.y_min is None and self.y_max is None:
            return array

        else:
            return array[self.z_min: self.z_max, self.y_min: self.y_max, self.x_min: self.x_max]

    def __copy__(self):
        return Box([self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max])

    @property
    def shape(self):
        if None in [self.x_max, self.x_min, self.y_max, self.y_min, self.z_min, self.z_max]:
            return None
        else:
            return self.z_min - self.z_max, self.y_max - self.y_min, self.x_max - self.x_min

    @classmethod
    def centered_at(cls, x0, y0, radius, z0=None):
        if z0 is None:
            return cls(x0 - radius, x0 + radius + 1, y0 - radius, y0 + radius + 1, None, None)
        else:
            return cls(x0 - radius, x0 + radius + 1, y0 - radius, y0 + radius + 1, z0 - radius, z0 + radius + 1)

    def crop_to_shape(self, shape):
        if self.x_min < 0:
            self.x_min = None
        if self.y_min < 0:
            self.y_min = None
        if self.z_min < 0:
            self.z_min = None
        if self.x_max >= shape[2]:
            self.x_max = None
        if self.y_max >= shape[1]:
            self.y_max = None
        if self.z_max >= shape[0]:
            self.z_max = None
