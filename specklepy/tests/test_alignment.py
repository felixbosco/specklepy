import numpy as np
import unittest

from specklepy.io.filearchive import FileArchive
from specklepy.core import alignment
from specklepy.plotting.utils import imshow


class TestAlignment(unittest.TestCase):

    def setUp(self):
        self.path = 'specklepy/tests/files/'
        self.files = FileArchive('synthetic/glao_600ms*.fits', in_dir=self.path).files
        self.shifts = [(0, 0), (34, -20), (-14, -51)]
        self.image_shape = (512, 512)
        self.cube_shape = (10, 512, 512)

    def test_get_pad_vectors(self):
        alignment.derive_pad_vectors(self.shifts)

    def test_pad_array(self):
        pad_vectors, ref_pad_vector = alignment.derive_pad_vectors(self.shifts)
        padded = alignment.pad_array(np.ones(self.image_shape), pad_vectors[1], mode='same',
                                     reference_image_pad_vector=ref_pad_vector)
        imshow(padded)
        pad_vectors, ref_pad_vector = alignment.derive_pad_vectors(self.shifts)
        for pad_vector in pad_vectors:
            padded = alignment.pad_array(np.ones(self.cube_shape), pad_vector=pad_vector, mode='same',
                                         reference_image_pad_vector=ref_pad_vector)


if __name__ == "__main__":
    unittest.main()
