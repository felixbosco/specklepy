from specklepy.utils.box import Box


class SubWindow(Box):

    @classmethod
    def from_str(cls, val):
        val.replace(' ', '')
        val.replace('[', '')
        val.replace(']', '')
        x_intv, y_intv = val.split(',')
        x_min, x_max = x_intv.split(':')
        y_min, y_max = y_intv.split(':')
        return cls(indexes=[x_min, x_max, y_min, y_max])
