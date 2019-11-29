import unittest
from holopy.io.filemanager import FileManager


class TestFileManager(unittest.TestCase):

    def setUp(self):
        self.filename = "data/test/example_cube.fits"

    def test_init(self):
        FileManager(self.filename)

    def test_str(self):
        fh = FileManager(self.filename)
        print(fh)

    def test_iter(self):
        fh = FileManager(self.filename)
        for file in fh:
            print(type(file), file)

if __name__ == "__main__":
    unittest.main()
