from astropy.units import Quantity


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

    def __init__(self, *args, scale, scaled=False):
        """

        Args:
            *args (any):
                Input for the x (and y) values of the tuple. Tuple and list-types are stored as x and y values. Int and
                float-types are stored and if only one value is parsed, then this is copied to both.
            scale (int or float or astropy.units.Quantity object):
                Scaling constant between pixel and the respective other physical unit.
            scaled (bool, optional):
                Indicate whether the positional arguments are provided in the scaled or unscaled space.
        """

        self.unit = None

        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (int, float)):
                x = arg
                y = arg
                unit = 1
            elif isinstance(arg, (list, tuple)):
                x = arg[0]
                y = arg[1]
            elif isinstance(arg, Quantity):
                if not scaled:
                    raise TypeError()

            else:
                raise TypeError(f"Type of positional argument not understood ({type(arg)})!")

        elif len(args) == 2:
            x = args[0]
            y = args[1]
        else:
            raise ValueError(f"ScaledTuple takes 1 or 2 positional arguments, not {len(args)}!")

        if scaled:
            self.x = x / scale
            self.y = y / scale
        else:
            self.x = x
            self.y = y

        self.scale = scale

    def __repr__(self):
        return f"ScaledTuple({self.x}, {self.y})"

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


class Position(ScaledTuple):

    def offset(self, other, scaled=False):
        if isinstance(other, tuple):
            if scaled:
                self.x += other[0] / self.scale
                self.y += other[1] / self.scale
            else:
                self.x += other[0]
                self.y += other[1]
        elif isinstance(other, ScaledTuple):
            self.x += other.x
            self.y += other.y


class ScaledShape(ScaledTuple):

    @property
    def shape(self):
        return self.index

    @property
    def scaled_shape(self):
        return self.coordinates

