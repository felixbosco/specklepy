from datetime import datetime
import glob
from IPython import embed
import numpy as np
import os
import string
import sys

from astropy.io import fits
from astropy.io.registry import IORegistryError
from astropy.table import Table

from specklepy.io import config
from specklepy.logging import logger
from specklepy.exceptions import SpecklepyTypeError


class FileArchive(object):

    """This class is handling the input files and is iterable.

    Attributes:
        file_list (str):
            Path to list of files or generic file path.
        ...
    """

    def __init__(self, file_list, in_dir=None, out_dir=None, **kwargs):
        """Create a FileArchive instance.

        Long description...

        Args:
            file_list (str, list):
                Path to list of files or generic file path. Can also be provided as list type.
            in_dir (str, optional):
                Path to the raw/ input data.
            out_dir (str, optional):
                Path to the product/ output data.
        """

        # Store in and out paths
        if in_dir is None:
            self.in_dir = './'
        else:
            self.in_dir = in_dir
        if out_dir is None:
            self.out_dir = './'
        else:
            self.out_dir = out_dir

        # Interpret the file list input
        if isinstance(file_list, str):
            # Search for files
            files = glob.glob(file_list)
            files.sort()
            if len(files) == 0:
                sys.tracebacklimit = 0
                raise FileNotFoundError("FileArchive did not find any file matching to{!r}.".format(file_list))
            else:
                logger.info("FileArchive found {} file(s) matching to {!r}.".format(len(files), file_list))

            if len(files) == 1 and not self.is_fits_file(files[0]):
                logger.info("Input file is not fits type. FileArchive assumes that input file {!r} contains file "
                            "names.".format(files[0]))
                self.table = self.read_table_file(files[0])

        elif isinstance(file_list, list):
            logger.info("FileArchive received a list of files.")
            self.table = self.gather_table_from_list(files=file_list, **kwargs)

        else:
            raise SpecklepyTypeError("FileArchive", 'file_list', type(file_list), 'str')

        # Log identified input files
        logger.debug("FileArchive lists the following files:")
        logger.debug(str(self.table))

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
        s = str(type(self)) + "\n"
        for file in self.files:
            s += "> {}\n".format(file)
        return s

    def __call__(self):
        return self.files

    def __len__(self):
        return len(self.files)

    @property
    def files(self):
        return self.table['FILE'].data

    @staticmethod
    def is_fits_file(filename):
        _, extension = os.path.splitext(filename)
        return extension == '.fits'

    def read_table_file(self, file, namekey='FILE'):
        """Interprets text in a file input."""

        try:
            table = Table.read(file)
        except IORegistryError:
            table = Table.read(file, format='ascii.fixed_width')
        except:
            files = []
            logger.info("Reading file names from input file {}.".format(file))
            with open(file, 'r') as f:
                for filename in f.readlines():
                    filename = filename.replace('\n', '')
                    files.append(filename)
            table = Table(data=[files], names=['FILE'], dtype=[object])

        return table

    @staticmethod
    def gather_table_from_list(files, cards, dtypes, names=None):
        """Gather file header information to fill the table

        Args:
            files (list):
                .
            cards (list):
                .
            dtypes (list):
                .
            names (list):
                .

        Returns:
            table (astropy.Table):
                .
        """

        # Initialize output file information table
        if names is None:
            table = Table(names=['FILE']+cards, dtype=[str]+dtypes)
        else:
            table = Table(names=['FILE']+names, dtype=[str]+dtypes)

        # Read data from files
        for file in files:
            logger.info(f"Retrieving header information from file {file}")
            hdr = fits.getheader(file)
            new_row = [os.path.basename(file)]
            for card in cards:
                try:
                    new_row.append(hdr[card])
                except KeyError:
                    logger.info(f"Skipping file {os.path.basename(file)} due to at least one missing header "
                                f"card ({card}).")
                    break
            if len(new_row) == len(table.columns):
                table.add_row(new_row)

        return table

    def write_table(self, file_name):
        logger.info(f"Writing file information to {file_name}")
        self.table.write(file_name, format='ascii.fixed_width', overwrite=True)

    def filter(self, filter_dict, namekey='FILE'):

        if not isinstance(filter_dict, dict):
            raise SpecklepyTypeError('FileArchive.filter', 'filter_dict', type(filter_dict), 'dict')

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

    def get_flats(self):
        return self.filter({'OBSTYPE': 'FLAT'})

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
        self.table['SETUP'] = [None] * len(self.table)

        # Iterate over keywords and identify unique settings per key
        for key in keywords:
            try:
                unique = np.unique(self.table[key].data)
            except KeyError:
                logger.info(f"Key {key} is not available in the file table and will be ignored!")
                continue
            logger.info(f"Identified {len(unique)} setups by keyword {key}:\t{unique}")

            for index, setup in enumerate(unique):
                for row in self.table:
                    if row[key] == setup:
                        if row['SETUP'] is None:
                            row['SETUP'] = str(index)
                        else:
                            row['SETUP'] + str(index)

        # Overwrite setup keys by length-1 string
        combinations = np.unique(self.table['SETUP'].data)
        for index, combination in enumerate(combinations):
            row_indexes = np.where(self.table['SETUP'].data == combination)
            self.table['SETUP'][row_indexes] = string.ascii_uppercase[index]

    def initialize_product_files(self, prefix='r'):
        product_files = []
        for file in self.filter({'OBSTYPE': 'SCIENCE'}):
            src = os.path.join(self.in_dir, file)
            dest = os.path.join(self.out_dir, prefix + file)
            logger.info(f"Initializing data product file {dest}")
            os.system(f"cp {src} {dest}")
            with fits.open(dest) as hdu_list:
                hdu_list[0].header.set('PIPELINE', 'Specklepy')
                hdu_list[0].header.set('REDUCED', datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))

            # Store new file in the list of product files
            product_files.append(dest)

        return product_files
