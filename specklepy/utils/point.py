class Point(object):

    def __init__(self, *args, order='xy'):

        self.x = None
        self.y = None
        self.z = None
        self.order = order

        for a, arg in enumerate(args):
            self.__setattr__(self.order[a], arg)

    def __repr__(self):
        return f"Point({self.x}, {self.y})"

    @property
    def ndim(self):
        return len(self.order)
