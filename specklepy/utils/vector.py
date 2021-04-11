import numpy as np


class Vector(object):

    def __init__(self, *args):
        if isinstance(args, tuple):
            self.values = list(args)
        if isinstance(args, list):
            self.values = args
        else:
            self.values = list(args)

    def __repr__(self):
        return f"Vector({self.values})"

    def __len__(self):
        return len(self.values)

    def __add__(self, other):
        if isinstance(other, (int, float)):
            from IPython import embed; embed()
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

    def __copy__(self):
        return Vector(*self.values)

    def length(self):
        return np.sqrt(self * self)
