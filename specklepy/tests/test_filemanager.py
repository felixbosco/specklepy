import unittest
from specklepy.io.filemanager import FileManager


class TestFileManager(unittest.TestCase):

    def setUp(self):
        self.filename = "specklepy/tests/files/example_cube.fits"
        self.filelist = "specklepy/tests/files/reduction/test_file_list.tab"

    def test_init(self):
        FileManager(self.filename)
        FileManager(self.filelist)

    def test_str(self):
        fh = FileManager(self.filename)
        print(fh)

    def test_iter(self):
        fh = FileManager(self.filename)
        for file in fh:
            print(type(file), file)

if __name__ == "__main__":
    unittest.main()
