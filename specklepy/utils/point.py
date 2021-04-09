import numpy as np


class Point(object):

    def __init__(self, *args, order='xy'):

        self.x = 0
        self.y = 0
        self.z = 0
        self.order = order

        for a, arg in enumerate(args):
            self.__setattr__(self.order[a], arg)

    def __repr__(self):
        return f"Point({self.x}, {self.y})"

    @property
    def ndim(self):
        return len(self.order)

    def __add__(self, other):
        dx = self.x + other.x
        dy = self.y + other.y
        dz = self.z + other.z
        return Point(dx, dy, dz, order='xyz')

    def __sub__(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return Point(dx, dy, dz, order='xyz')

    def __len__(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    def distance_to(self, other):
        return len(self - other)
