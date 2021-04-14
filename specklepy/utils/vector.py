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
            raise TypeError(f"unsupported operand type(s) for +: 'Vector' and {type(other)!r}")

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.subtract(self.values, other))
        elif isinstance(other, Vector):
            return Vector(np.subtract(self.values, other.values))
        else:
            raise TypeError(f"unsupported operand type(s) for -: 'Vector' and {type(other)!r}")

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.multiply(self.values, other))
        elif isinstance(other, Vector):
            return np.sum(np.multiply(self.values, other.values))
        else:
            raise TypeError(f"unsupported operand type(s) for *: 'Vector' and {type(other)!r}")

    def __floordiv__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.floor_divide(self.values, other))
        else:
            raise TypeError(f"unsupported operand type(s) for //: 'Vector' and {type(other)!r}")

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return Vector(np.true_divide(self.values, other))
        else:
            raise TypeError(f"unsupported operand type(s) for /: 'Vector' and {type(other)!r}")

    def __getitem__(self, item):
        return self.values.__getitem__(item)

    def __setitem__(self, item, value):
        return self.values.__setitem__(item, value)

    def __copy__(self):
        return Vector(*self.values)

    def length(self):
        return np.sqrt(self * self)

    @property
    def dtype(self):
        return type(self.values[0])

    @dtype.setter
    def dtype(self, type):
        self.values = [type(x) for x in self.values]

    def astype(self, type):
        return type(self.values)

    def round(self):
        self.values = [round(x) for x in self.values]
