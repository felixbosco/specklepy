from IPython import embed
import numpy as np

from specklepy.utils.box import Box


class SubWindow(Box):

    @classmethod
    def from_str(cls, val=None):

        # Create full box window if no string provided
        if val is None:
            return cls(indexes=None)

        # Unravel coordinates from string
        val = val.replace(' ', '')
        val = val.replace('[', '')
        val = val.replace(']', '')
        x_intv, y_intv = val.split(',')
        x_min, x_max = x_intv.split(':')
        y_min, y_max = y_intv.split(':')
        return cls(indexes=np.array([x_min, x_max, y_min, y_max]).astype(int))
