class SpecklepyValueError(ValueError):

    """ValueError with message that is composed in a default format."""

    def __init__(self, objname, argname, argvalue, expected, *args, **kwargs):
        """Raise a ValueError with formatted message.

        Args:
            objname (str):
                Name of the object/ function/ etc. that is raising the error.
            argname (str):
                Name of the argument that received input of value.
            argvalue (str):
                Value that was received.
            expected (str):
                Expected value for argument argname.
            *args:
                Passed to ValueError.
            **kwargs:
                Passed to ValueError.
        """
        message = "{!r} received argument {}={}, but must be {!r}!".format(objname, argname, argvalue, expected)
        # super().__init__(message, *args, **kwargs)
        raise ValueError(message, *args, **kwargs)


class SpecklepyTypeError(TypeError):

    """TypeError with message that is composed in a default format."""

    def __init__(self, objname, argname, argtype, expected, *args, **kwargs):
        """Raise a TypeError with formatted message.

        Args:
            objname (str):
                Name of the object/ function/ etc. that is raising the error.
            argname (str):
                Name of the argument that received input of value.
            argtype (str, type, any):
                Type that was received by object. If not str or type, the type will be evaluated.
            expected (str):
                Expected type for argument argname.
            *args:
                Passed to TypeError.
            **kwargs:
                Passed to TypeError.
        """

        if not isinstance(argtype, (str, type)):
            argtype = type(argtype)

        message = "{!r} received argument {!r} of type {!r}, but must be {!r} type!".format(objname, argname, argtype, expected)
        # super().__init__(message, *args, **kwargs)
        raise TypeError(message, *args, **kwargs)
