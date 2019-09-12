import matplotlib.pyplot as plt
from holopy.config import plotting

def imshow(image, title=None, norm=None):
    plt.figure()
    plt.imshow(image, norm=norm)
    plt.title(title)
    plt.colorbar(pad=0.0)
    plt.show()
    plt.close()

def plot_powerspec1d(image):
    pass
