class SpecklepyValueError(ValueError):

    def __init__(self, object, argname, argvalue, expected, *args, **kwargs):
        super().__init__(*args, **kwargs)
        message = "{!r} received argument {}={}, but must be {!r}!".format(object, argname, argvalue, expected)
        raise ValueError(message, *args, **kwargs)


class SpecklepyTypeError(TypeError):

    def __init__(self, object, argname, argtype, expected, *args, **kwargs):
        super().__init__(*args, **kwargs)
        message = "{!r} received argument {!r} of type {!r}, but must be {!r} type!".format(object, argname, argtype, expected)
        raise TypeError(message, *args, **kwargs)
