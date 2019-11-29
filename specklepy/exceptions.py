class SpecklepyValueError(Exception):

    def __init__(self, function, argument, value, expected):
        message = "The function {} received {}={}, but must be {}!".format(function, argument, value, expected)
        raise ValueError(message)


class SpecklepyTypeError(Exception):

    def __init__(self, function, argument, type, expected):
        message = "The function {} received argument {} of type {}, but must be {}".format(function, argument, value, expected)
        raise TypeError(message)
