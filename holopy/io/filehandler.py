import glob
from astropy.io import fits

from holopy.logging import logging
from holopy.io.outfile import Outfile


class FileHandler(object):

    """
    This class is handling the input files and is iterable.
    """

    def __init__(self, input):
        self._characterize_input(input)
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

    def __str__(self):
        s = "<holopy.io.filehandler.FileHandler object>\n"
        s += str(self.files)
        return s


    def _characterize_input(self, input):
        """
        This function interpretes the file input and stores a list to self.files.
        """
        if isinstance(input, str):
            if '*' in input:
                logging.info("Input is a generic file name.")
                self.input_type = 'generic'
                # find generic files with glob.glob()
                self.files = glob.glob(input)
                logging.info("Found {} files matching to the generic name {}.".format(len(self.files), input))
            else:
                extension = input.split('.')[-1]
                if extension == 'fits':
                    try:
                        fits.getheader(input)
                    except FileNotFoundError as e:
                        logging.warning("Fits file <{}> has not been found.".format(input))
                    self.files = [input]
                elif extension == 'spam':
                    logging.info("Input is the name of a spam spectrum file.")
                    self.files = [input]
                else:
                    logging.info("Assuming that input file {} contains file names.".format(input))
                    self._read_file_list_file(input)
        elif isinstance(input, list):
            logging.info("Input is a list of files.")
            self.input_type = 'list'
            self.files = input
        else:
            raise TypeError("FileHandler got input of unexpected type ({}).".format(type(input)))


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


    def create_outfile(self):
        header = fits.getheader(self.files[0])
        return Outfile(hdr=header, file_list=self.files)
