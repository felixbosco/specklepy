import unittest
from specklepy.io.filearchive import FileArchive
from specklepy.deprecated.reconstructionfile import ReconstructionFile


class TestReconstructionFile(unittest.TestCase):

    def setUp(self):
        self.FileManager = FileArchive("specklepy/tests/files/example_cube.fits")

    def test_init(self):
        ReconstructionFile(filename="specklepy/tests/files/test_reconstructionfile.fits", files=self.FileManager.files, cards={"RECONSTRUCTION": "Test"})

    def test_set_data(self):
        pass

if __name__ == "__main__":
    unittest.main()
