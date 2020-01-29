import os
import glob
from astropy.io.registry import IORegistryError
from astropy.table import Table

from specklepy.logging import logging


class FileManager(object):

    """
    This class is handling the input files and is iterable.
    """

    def __init__(self, input):
        """Instantiate a FileManager object.

        Long description...

        Args:
            input (str):

        """

        self.input = input
        if isinstance(self.input, str):
            # Search for files
            self.files = glob.glob(self.input)
            self.files.sort()
            logging.info("FileManager found {} file(s) matching to {!r}.".format(len(self.files), input))

            if len(self.files) == 1 and not self.is_fits_file(self.files[0]):
                logging.info("Input file is not fits type. FileManager assumes that input file {!r} contains file names.".format(self.files[0]))
                self.extract_file_names(self.files[0])

        elif isinstance(input, list):
            logging.info("FileManager received a list of files.")
            self.files = input

        else:
            raise TypeError("FileManager received input of unexpected type ({}).".format(type(input)))

        # Log identified input files
        logging.info("FileManager lists the following files:")
        for f, file in enumerate(self.files):
            logging.info("{:4d}: {}".format(f+1, file))

        # Initialize the index for iteration
        self.index = 0


    def __iter__(self):
        return self

    def __next__(self):
        try:
            result = self.files[self.index]
        except IndexError:
            raise StopIteration
        self.index += 1
        return result

    def __getitem__(self, index):
        return self.files[index]

    def __setitem__(self, index, value):
        self.files[index] = value

    def __str__(self):
        s = "<specklepy.io.filemanager.FileManager object>\n"
        s += "> List of files:\n"
        for file in self.files:
            s += "> {}\n".format(file)
        return s

    def __call__(self):
        return self.files

    def __len__(self):
        return len(self.files)


    def is_fits_file(self, filename):
        _, extension = os.path.splitext(filename)
        return extension == '.fits'


    def extract_file_names(self, file, namekey='FILE'):
        """
        Interpretes text in a file input.
        """

        try:
            self.filelist = Table.read(file)
            self.files = self.filelist[namekey].data
        except IORegistryError:
            self.filelist = Table.read(file, format='ascii.fixed_width')
            self.files = self.filelist[namekey].data
        except:
            self.files = []
            logging.info("Reading file names from input file {}.".format(file))
            with open(file, 'r') as f:
                for filename in f.readlines():
                    filename = filename.replace('\n', '')
                    self.files.append(filename)
