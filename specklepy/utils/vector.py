import numpy as np


class Vector(object):

    def __init__(self, *args):
        args = list(args)
        if len(args) == 1:
            if isinstance(args[0], (int, float)):
                self.values = args
            if isinstance(args[0], (tuple, list, np.ndarray)):
                self.values = list(args[0])
            else:
                raise TypeError
        else:
            self.values = args

    def __repr__(self):
        return f"Vector({self.values})"

    def __len__(self):
        return len(self.values)

    def __add__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.add(self.values, other))
        elif isinstance(other, Vector):
            return Vector(np.add(self.values, other.values))
        else:
            raise TypeError

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.subtract(self.values, other))
        elif isinstance(other, Vector):
            return Vector(np.subtract(self.values, other.values))
        else:
            raise TypeError

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.multiply(self.values, other))
        elif isinstance(other, Vector):
            return np.sum(np.multiply(self.values, other.values))
        else:
            raise TypeError

    def __floordiv__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.floor_divide(self.values, other))
        else:
            raise TypeError

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.true_divide(self.values, other))
        else:
            raise TypeError

    def __copy__(self):
        return Vector(*self.values)

    def length(self):
        return np.sqrt(self * self)
