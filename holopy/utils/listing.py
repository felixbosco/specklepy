class Listing(object):

    def __init__(self, **kwargs):
        for key in kwargs:
            self.__setattr__(key, kwargs[key])

    def __str__(self):
        return str(self.__dict__)
