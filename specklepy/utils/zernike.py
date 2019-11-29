import numpy as np
from math import factorial
import warnings


class Zernike(object):

    index_schemes = {
        'Noll': [(0,0), (1,1), (1,-1), (2,0), (2,-2), (2,2), (3,-1), (3,1), (3,-3), (3,3), (4,0), (4,2), (4,-2), (4,4), (4,-4), (5,1), (5,-1), (5,3), (5,-3), (5,5)],
        'OSA': [(0,0), (1,-1), (1,1), (2,-2), (2,0), (2,2), (3,-3), (3,-1), (3,1), (3,3), (4,-4), (4,-2), (4,0), (4,2), (4,4), (5,-5), (5,-3), (5,-1), (5,1), (5,3)],
        'Fringe': [(0,0), (1,1), (1,-1), (2,0), (2,2), (2,-2), (3,1), (3,-1), (4,0), (3,3), (3,-3), (4,2), (4,-2), (5,1), (5,-1), (6,0), (4,4), (4,-4), (5,3), (5,-3)]
    }

    classical_names = {
        'piston': (0,0),
        'tilt': (1,-1),
        'tip': (1,1),
        'defocus': (2,0)
    }

    def __init__(self, index_scheme='Noll'):
        self.index_scheme = index_scheme
        if index_scheme == 'ANSI':
            index_scheme = 'OSA'
        self.indices = self.index_schemes[index_scheme]

    def __call__(self, coeffs, size=32):
        rho = self.init_rho(size=size)
        phi = self.init_phi(size=size)
        out = np.zeros((size, size))

        # Iterate over indices according to index scheme
        for index, coeff in enumerate(coeffs):
            n, m = self.indices[index]
            out += coeff * self.evaluate(n, m, rho=rho, phi=phi)

        mask = np.ma.masked_greater(rho, size / 2).mask
        return np.ma.masked_array(out, mask=mask)

    def init_rho(self, size):
        xx, yy = np.mgrid[:size, :size] - (size / 2)
        rho = np.sqrt(np.square(xx) + np.square(yy))
        rho = rho / (size / 2) # Normalize to order unity
        return np.ma.masked_greater_equal(rho, 1.0)

    def init_phi(self, size):
        xx, yy = np.mgrid[:size, :size] - (size / 2)
        return np.angle(xx*1j + yy)

    def evaluate(self, n, m, rho, phi):
        if n < abs(m):
            raise ValueError('n must be greater or equal to m!')
        if m >= 0 :
            return self.radial(n, abs(m), rho) * np.cos(m * phi)
        else:
            return self.radial(n, abs(m), rho) * np.sin(m * phi)

    def radial(self, n, m, rho):
        if (n - m) % 2: # odd
            return np.zeros(rho.shape)
        else: # even
            kmax = int((n - m) / 2)
            out = np.zeros(rho.shape)
            for k in range(kmax + 1):
                factorials = np.divide(np.power(-1, k) * factorial(n - k), factorial(k) * factorial((n+m)/2 - k) * factorial((n-m)/2 - k))
                out += factorials * np.power(rho, n - 2 * k)
            return np.ma.masked_array(out, mask=rho.mask)

    #
    def defocus(self, coeff, size):
        return coeff * self.evaluate(*self.classical_names['defocus'], self.init_rho(size), self.init_phi(size))
