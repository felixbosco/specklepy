import numpy as np

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.logging import logger


class Box(object):

    def __init__(self, x_min=None, x_max=None, y_min=None, y_max=None, indexes=None):

        # Initialize attributes with kwargs
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        # Interpret indexes input and overwrite
        if indexes is not None:
            if not isinstance(indexes, (tuple, list)):
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

    def __repr__(self):
        return f"[[{self.y_min}, {self.y_max}], [{self.x_min}, {self.x_max}]]"

    def get(self, array):

        # Abort if Box is undefined
        if self.x_min is None and self.x_max is None and self.y_min is None and self.y_max is None:
            return array

        else:
            # Set limits
            x_min = 0 if self.x_min is None else self.x_min
            x_max = array.shape[1] if self.x_max is None else self.x_max
            y_min = 0 if self.y_min is None else self.y_min
            y_max = array.shape[1] if self.y_max is None else self.y_max

            # Create index sets
            x, y = np.meshgrid(range(y_min, y_max), range(x_min, x_max))

            return array[y, x]
