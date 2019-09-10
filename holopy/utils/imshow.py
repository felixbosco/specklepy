import matplotlib.pyplot as plt

def imshow(image, norm=None):
    plt.imshow(image, norm=norm)
    plt.colorbar(pad=0.0)
    plt.show()
    plt.close()
