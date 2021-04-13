import numpy as np
import unittest

from specklepy.io.filearchive import FileArchive
from specklepy.core.alignment import FrameAlignment
from specklepy.plotting.utils import imshow


class TestAlignment(unittest.TestCase):

    def setUp(self):
        self.path = 'specklepy/tests/files/'
        self.files = FileArchive('synthetic/glao_600ms*.fits', file_path=self.path).files
        self.shifts = [(0, 0), (34, -20), (-14, -51)]
        self.image_shape = (512, 512)
        self.cube_shape = (10, 512, 512)

    def test_get_pad_vectors(self):
        alignment = FrameAlignment()
        alignment.derive_pad_vectors(self.shifts)

    def test_pad_array(self):
        alignment = FrameAlignment()
        pad_vectors, ref_pad_vector = alignment.derive_pad_vectors(self.shifts)
        padded = alignment.pad_array(np.ones(self.image_shape), pad_vector_index=1, mode='same')
        imshow(padded)
        pad_vectors, ref_pad_vector = alignment.derive_pad_vectors(self.shifts)
        for p, pad_vector in enumerate(pad_vectors):
            padded = alignment.pad_array(np.ones(self.cube_shape), pad_vector_index=p, mode='same')


if __name__ == "__main__":
    unittest.main()
