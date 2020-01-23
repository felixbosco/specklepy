import os
import glob
from astropy.io import fits

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
        if isinstance(input, str):
            self.files = glob.glob(input)
            self.files.sort()
            logging.info("FileManager found {} file(s) matching to {}.".format(len(self.files), input))

            for file in self.files:
                extension = file.split('.')[-1]
                if extension == 'fits':
                    try:
                        fits.getheader(file)
                    except FileNotFoundError as e:
                        logging.warning("Fits file <{}> has not been found.".format(file))
                    # self.files = [input]
                else:
                    logging.info("Input file is not fits type. FileManager assumes that input file {} contains file names.".format(input))
                    self._read_file_list_file(input)
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


    def _characterize_input(self, input):
        """
        This function interpretes the file input and stores a list to self.files.
        """
        if isinstance(input, str):
            if '*' in input or '?'in input:
                logging.info("Input is a generic file name.")
                self.input_type = 'generic'
                # find generic files with glob.glob()
                self.files = glob.glob(input)
                self.files.sort()
                logging.info("Found {} files matching to the generic name {}.".format(len(self.files), input))
            else:
                extension = input.split('.')[-1]
                if extension == 'fits':
                    try:
                        fits.getheader(input)
                    except FileNotFoundError as e:
                        logging.warning("Fits file <{}> has not been found.".format(input))
                    self.files = [input]
                # elif extension == 'spam':
                #     logging.info("Input is the name of a spam spectrum file.")
                #     self.files = [input]
                else:
                    logging.info("Assuming that input file {} contains file names.".format(input))
                    self._read_file_list_file(input)
        elif isinstance(input, list):
            logging.info("Input is a list of files.")
            self.input_type = 'list'
            self.files = input
        else:
            raise TypeError("FileManager got input of unexpected type ({}).".format(type(input)))


    def _read_file_list_file(self, input):
        """
        Interpretes text file input.
        """
        self.files = []
        logging.info("Reading files from input file list {}.".format(input))
        with open(input, 'r') as f:
            for filename in f.readlines():
                filename = filename.replace('\n', '')
                self.files.append(filename)
