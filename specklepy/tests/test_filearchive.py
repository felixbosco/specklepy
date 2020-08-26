import unittest
from specklepy.io.filearchive import FileArchive


class TestFileArchive(unittest.TestCase):

    def setUp(self):
        self.filename = "specklepy/tests/files/example_cube.fits"
        self.filelist = "specklepy/tests/files/reduction/test_file_list.tab"

    def test_init(self):
        FileArchive(self.filename)
        FileArchive(self.filelist)

    def test_str(self):
        fh = FileArchive(self.filename)
        print(fh)

    def test_iter(self):
        fh = FileArchive(self.filename)
        for file in fh:
            print(type(file), file)

if __name__ == "__main__":
    unittest.main()
