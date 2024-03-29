import numpy as np

from specklepy.utils.box import Box


class SubWindow(Box):

    @classmethod
    def from_str(cls, s=None, full=None, order='yx'):

        # Create full box window if no string provided
        if s is None:
            return cls(indexes=None)

        # Unravel coordinates from string
        indexes = cls.unravel_string(s=s, order=order)

        # Provide indexes relative to a full window
        if full is not None:
            full = cls.unravel_string(s=full, order=order)
            indexes = indexes - full
            indexes = np.where(indexes == 0, None, indexes)  # Avoid IndexErrors by substituting the edges by Nones

        return cls(indexes=indexes)

    @staticmethod
    def unravel_string(s, order='yx'):
        s = s.replace(' ', '')
        s = s.replace('[', '')
        s = s.replace(']', '')
        if order == 'xy':
            x_intv, y_intv = s.split(',')
        elif order == 'yx':
            y_intv, x_intv = s.split(',')
        else:
            raise ValueError(f"Order {order!r} not understood in unraveling sub-window indexes!")
        x_min, x_max = x_intv.split(':')
        y_min, y_max = y_intv.split(':')
        return np.array([x_min, x_max, y_min, y_max]).astype(int)

    # def switch_to_zero_based_indexing(self):
    #     self.x_min -= 1
    #     # self.x_max -= 1
    #     self.y_min -= 1
    #     # self.y_max -= 1
