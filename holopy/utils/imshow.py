import matplotlib.pyplot as plt

def imshow(image):
    plt.imshow(image)
    plt.colorbar(pad=0.0)
    plt.show()
    plt.close()
