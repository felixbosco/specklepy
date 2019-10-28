import numpy as np
from zernike import Zernike
from plot import imshow

z = Zernike()

if False:
    rho = z.init_rho(256)
    imshow(rho, title='Radius')

if False:
    phi = z.init_phi(256)
    imshow(phi, title='Azimuth')

if True:
    coeffs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    # coeffs = np.random.rand((10))
    out = z(coeffs, size=128)
    imshow(out, title='Zernike polynomial {}'.format(coeffs))

if False:
    imshow(z.defocus(-1, 256), title='defocus')
