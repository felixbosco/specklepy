import numpy as np

def get_max(array):
    return np.unravel_index(np.argmax(array, axis=None), array.shape)
