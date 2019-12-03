import numpy as np
from astropy.io import fits



def subtract(params, debug=False):
    print(params.inFiles)
    print(params.skyFiles)
