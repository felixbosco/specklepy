import numpy as np
import matplotlib.pyplot as plt

class MatchStarList(object):

    def __init__(self):
        pass

    def __call__(self, files):
        for index, star_list in enumerate(files):
            table = np.loadtxt(star_list, skiprows=1).transpose()
            plt.plot(table[0], table[1], '.')
        plt.show()
        plt.close()
