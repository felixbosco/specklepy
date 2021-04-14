from .box import Box, Box3D
from .point import Point
from .time import default_time_stamp
from .vector import Vector


def save_eval(val):
    """Evaluate variables, assuming they are str representations of e.g. int-types.

    Args:
        val (any):
            Value that shall be evaluated.

    Returns:
        res (any):
            Evaluated version of value if built-in `eval` was successful, otherwise `val`.
    """

    try:
        res = eval(val)
    except TypeError:
        # val is not str-type, so return
        res = val
    except NameError:
        # val is str-type, but not the name of a variable, so return
        res = val
    except SyntaxError:
        # val is str-type, but does not just contain a representation, so return
        res = val

    return res
