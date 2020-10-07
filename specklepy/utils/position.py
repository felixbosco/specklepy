class Position(object):

    """A class for transforming positions between pixel coordinates and unit coordinates.

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
            self.center = (0, 0)
        else:
            self.center = center

    @property
    def x(self):
        return self._x - self.center[0]

    @property
    def y(self):
        return self._y - self.center[1]

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
