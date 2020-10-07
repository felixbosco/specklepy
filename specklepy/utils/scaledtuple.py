import numpy as np


class ScaledTuple(object):

    """A class for transforming tuples (positions and shapes) between pixel coordinates and unit coordinates.

    Attributes:
         x (float)
         y (float)
         x_scaled (float)
         y_scaled (float)
         index (tuple)
         coordinates (tuple)
         center (tuple)
    """

    def __init__(self, x, y, scale, scaled=False, center=None):

        if scaled:
            self._x = x / scale
            self._y = y / scale
        else:
            self._x = x
            self._y = y

        self.scale = scale

        if center is None:
            self._center = (0, 0)
        else:
            if scaled:
                self._center = tuple(np.array(center) / scale)
            else:
                self._center = center

    def __repr__(self):
        return f"ScaledTuple({self._x}, {self._y})"

    @property
    def x(self):
        return self._x - self._center[0]

    @property
    def y(self):
        return self._y - self._center[1]

    @property
    def index(self):
        return int(self.x), int(self.y)

    @property
    def x_scaled(self):
        return self.x * self.scale

    @property
    def y_scaled(self):
        return self.y * self.scale

    @property
    def coordinates(self):
        return self.x_scaled, self.y_scaled

    def shift(self, other, scaled=False):
        if isinstance(other, tuple):
            if scaled:
                self._x += other[0] / self.scale
                self._y += other[1] / self.scale
            else:
                self._x += other[0]
                self._y += other[1]
        elif isinstance(other, Position):
            self._x += other._x
            self._y += other._y
