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

    def __repr__(self):
        return f"[[{self.x_min}, {self.x_max}], [{self.y_min}, {self.y_max}]]"

    def __call__(self, array):  #, update=False):

        # Abort if Box is undefined
        if self.x_min is None and self.x_max is None and self.y_min is None and self.y_max is None:

            # # Update limits
            # if update:
            #     self.x_min = 0
            #     self.x_max = array.shape[1] - 1
            #     self.y_min = 0
            #     self.y_max = array.shape[0] - 1

            return array

        else:
            # # Set limits
            # x_min = 0 if self.x_min is None else self.x_min
            # x_max = array.shape[1] - 1 if self.x_max is None else self.x_max
            # y_min = 0 if self.y_min is None else self.y_min
            # y_max = array.shape[0] - 1 if self.y_max is None else self.y_max
            #
            # # Update limits
            # if update:
            #     self.x_min = x_min
            #     self.x_max = x_max
            #     self.y_min = y_min
            #     self.y_max = y_max
            #
            # # Create index sets
            # x, y = np.meshgrid(range(self.y_min, self.y_max), range(self.x_min, self.x_max))

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

        return obj

    def crop_to_shape(self, shape):
        if self.x_min < 0:
            self.x_min = 0
        if self.y_min < 0:
            self.y_min = 0
        if self.x_max >= shape[1]:
            self.x_max = -1
        if self.y_min >= shape[0]:
            self.y_min = -1

    def transpose(self):
        x_min = self.x_min
        x_max = self.x_max
        self.x_min = self.y_min
        self.x_max = self.y_max
        self.y_min = x_min
        self.y_max = x_max
