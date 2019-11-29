import unittest
import numpy as np
from specklepy.io.filemanager import FileManager
from specklepy.core import alignment
from specklepy.utils.plot import imshow


class Testalignment(unittest.TestCase):

    def setUp(self):
        self.files = FileManager('/home/bosco/Documents/sowat/simulations/noAO_200ms_x100*_st4_shift.fits')()
        self.shifts = [(0, 0), (34, -20), (-14, -51)]
        self.image_shape = (1024, 1024)
        self.cube_shape = (100, 1024, 1024)

    def test_get_shifts(self):
        alignment.get_shifts(self.files)

    def test_get_pad_vectors(self):
        # Test all modes
        for mode in ['same', 'full', 'valid']:
            alignment.get_pad_vectors(self.shifts, self.image_shape, self.image_shape, mode=mode)
        # Test switching to cube mode
        alignment.get_pad_vectors(self.shifts, self.cube_shape, self.image_shape, mode='same')
        # Test error handling
        with self.assertRaises(ValueError):
            alignment.get_pad_vectors(self.shifts, self.image_shape, self.image_shape, mode='Nonsense')

    def test_pad_array(self):
        pad_vectors, ref_pad_vector = alignment.get_pad_vectors(self.shifts, self.image_shape, self.image_shape, mode='same')
        padded = alignment.pad_array(np.ones(self.image_shape), pad_vectors[1], mode='same', reference_image_pad_vector=ref_pad_vector)
        imshow(padded)
        pad_vectors, ref_pad_vector = alignment.get_pad_vectors(self.shifts, self.cube_shape, self.image_shape, mode='same')
        for pad_vector in pad_vectors:
            padded = alignment.pad_array(np.ones(self.cube_shape), pad_vector=pad_vector, mode='same', reference_image_pad_vector=ref_pad_vector)
            print('>>>', padded.shape)



if __name__ == "__main__":
    unittest.main()
