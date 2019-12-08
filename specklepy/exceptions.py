class SpecklepyValueError(ValueError):

    def __init__(self, object, argname, argvalue, expected, *args, **kwargs):
        super().__init__(*args, **kwargs)
        message = "{} received argument {}={}, but must be {}!".format(object, argname, argvalue, expected)
        raise ValueError(message, *args, **kwargs)


class SpecklepyTypeError(TypeError):

    def __init__(self, object, argname, argtype, expected, *args, **kwargs):
        super().__init__(*args, **kwargs)
        message = "{} received argument {} of type {}, but must be {} type!".format(object, argname, argtype, expected)
        raise TypeError(message, *args, **kwargs)
