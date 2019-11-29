import unittest
from holopy.io.filemanager import FileManager
from holopy.core import alignment


class Testalignment(unittest.TestCase):

    def setUp(self):
        self.files = FileManager('/home/bosco/Documents/sowat/simulations/noAO_200ms_x100*_st4_shift.fits')()
        self.shifts = [(0, 0), (34, -20), (-14, -51)]
        self.image_shape = (1024, 1024)


    def test_compute_shifts(self):
        alignment.compute_shifts(self.files)

    def test_get_pad_vectors(self):
        alignment.get_pad_vectors(self.shifts, self.image_shape, self.image_shape, mode='same')
        alignment.get_pad_vectors(self.shifts, self.image_shape, self.image_shape, mode='full')
        alignment.get_pad_vectors(self.shifts, self.image_shape, self.image_shape, mode='valid')
        with self.assertRaises(ValueError):
            alignment.get_pad_vectors(self.shifts, self.image_shape, self.image_shape, mode='Nonsense')


if __name__ == "__main__":
    unittest.main()
