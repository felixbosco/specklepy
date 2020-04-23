import glob
import os
import string
import numpy as np
from astropy.io.registry import IORegistryError
from astropy.table import Table

from specklepy.logging import logger
from specklepy.exceptions import SpecklepyTypeError


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
            if len(self.files) == 0:
                raise FileNotFoundError("FileManager did not find any file matching to{!r}.".format(input))
            else:
                logger.info("FileManager found {} file(s) matching to {!r}.".format(len(self.files), input))

            if len(self.files) == 1 and not self.is_fits_file(self.files[0]):
                logger.info("Input file is not fits type. FileManager assumes that input file {!r} contains file names.".format(self.files[0]))
                self.extract_file_names(self.files[0])

        elif isinstance(input, list):
            logger.info("FileManager received a list of files.")
            self.files = input

        else:
            raise SpecklepyTypeError("FileManager", 'input', type(input), 'str')# received input of unexpected type ({}).".format(type(input)))

        # Log identified input files
        logger.debug("FileManager lists the following files:")
        for f, file in enumerate(self.files):
            logger.debug("{:4d}: {}".format(f+1, file))

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
            self.table = Table.read(file)
            self.files = self.table[namekey].data
        except IORegistryError:
            self.table = Table.read(file, format='ascii.fixed_width')
            self.files = self.table[namekey].data
        except:
            self.files = []
            logger.info("Reading file names from input file {}.".format(file))
            with open(file, 'r') as f:
                for filename in f.readlines():
                    filename = filename.replace('\n', '')
                    self.files.append(filename)

        # Replace finite length strings of file names by object type for
        # file names of arbitrary lengths
        file_column = np.array(self.table[namekey], dtype=object)
        self.table[namekey] = file_column


    def filter(self, filter_dict, namekey='FILE'):

        if not isinstance(filter_dict, dict):
            raise SpecklepyTypeError('FileManager.filter', 'filter_dict', type(filter_dict), 'dict')

        mask = [True] * len(self.table[namekey])
        for index, key in enumerate(filter_dict.keys()):
            if isinstance(filter_dict[key], list):
                submask = [False] * len(self.table[namekey])
                for correct in filter_dict[key]:
                    submask |= (self.table[key] == correct)
                mask &= submask
            else:
                mask &= (filter_dict[key] == self.table[key])

        return self.table[namekey][mask].data


    def update_filenames(self, filenames):
        """Updates the file names in the files table.

        Args:
            filenames (dict):
                Dictionary mapping the current file names onto new ones.
        """

        if not isinstance(filenames, dict):
            raise SpecklepyTypeError('FileManager.update_filenames', argname='filenames', argtype=type(filenames), expected='dict')

        logger.info('Updating file names:')
        for outdated in filenames.keys():
            index = np.where(self.table['FILE'] == outdated)
            self.table['FILE'][index] = filenames[outdated]
            logger.info(f"{outdated} > {filenames[outdated]}")


    def identify_setups(self, keywords):
        """Identify distinct observational setups in a list of files.

        This function identifies distinct observational setups.

        Args:
            keywords (list of str):
        """

        # Check input parameters
        if not isinstance(keywords, list):
            raise SpecklepyTypeError('identify_setups', 'keywords', type(keywords), 'list')

        # Identifying setups key-by-key
        logger.info("Identifying distinct observational setups in the file list...")
        self.table['Setup'] = [None] * len(self.table['FILE'])

        for key in keywords:
            try:
                unique = np.unique(self.table[key].data)
            except KeyError:
                logger.info(f"Key {key} is not available in the file table and will be ignored!")
                continue
            logger.info(f"Identified {len(unique)} setups by keyword {key}:")
            logger.info(f"\t{unique}")
            if len(unique) == 1:
                continue

            for index, setup in enumerate(unique):
                for row in self.table:
                    if row[key] == setup:
                        if row['Setup'] is None:
                            row['Setup'] = [str(index)]
                        else:
                            row['Setup'].append(str(index))

        for row in self.table:
            row['Setup'] = ''.join(row['Setup'])

        # Overwrite setup keys by length-1 string
        combinations = np.unique(self.table['Setup'].data)
        for index, combination in enumerate(combinations):
            row_indizes = np.where(self.table['Setup'].data == combination)
            self.table['Setup'][row_indizes] = string.ascii_uppercase[index]
