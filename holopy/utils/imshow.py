import matplotlib.pyplot as plt

def imshow(image, title=None, norm=None):
    plt.figure()
    plt.imshow(image, norm=norm)
    plt.title(title)
    plt.colorbar(pad=0.0)
    plt.show()
    plt.close()
