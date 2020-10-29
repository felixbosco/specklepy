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
            if not isinstance(indexes, list):
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

    def __call__(self, array, update=True):

        # Abort if Box is undefined
        if self.x_min is None and self.x_max is None and self.y_min is None and self.y_max is None:

            # Update limits
            if update:
                self.x_min = 0
                self.x_max = array.shape[1] - 1
                self.y_min = 0
                self.y_max = array.shape[0] - 1

            return array

        else:
            # Set limits
            x_min = 0 if self.x_min is None else self.x_min
            x_max = array.shape[1] - 1 if self.x_max is None else self.x_max
            y_min = 0 if self.y_min is None else self.y_min
            y_max = array.shape[0] - 1 if self.y_max is None else self.y_max

            # Update limits
            if update:
                self.x_min = x_min
                self.x_max = x_max
                self.y_min = y_min
                self.y_max = y_max

            # Create index sets
            x, y = np.meshgrid(range(y_min, y_max), range(x_min, x_max))

            return array[y, x]

    @property
    def shape(self):
        if None in [self.x_max, self.x_min, self.y_max, self.y_min]:
            return None
        else:
            return self.x_max - self.x_min, self.y_max - self.y_min
