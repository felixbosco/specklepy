from datetime import datetime
import glob
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

    def __init__(self, file_list, in_dir=None, out_dir=None):
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
            self.files = glob.glob(file_list)
            self.files.sort()
            if len(self.files) == 0:
                sys.tracebacklimit = 0
                raise FileNotFoundError("FileArchive did not find any file matching to{!r}.".format(file_list))
            else:
                logger.info("FileArchive found {} file(s) matching to {!r}.".format(len(self.files), file_list))

            if len(self.files) == 1 and not self.is_fits_file(self.files[0]):
                logger.info("Input file is not fits type. FileArchive assumes that input file {!r} contains file "
                            "names.".format(self.files[0]))
                self.read_table_file(self.files[0])

        elif isinstance(file_list, list):
            logger.info("FileArchive received a list of files.")
            self.files = file_list

        else:
            raise SpecklepyTypeError("FileArchive", 'file_list', type(file_list), 'str')

        # Log identified input files
        logger.debug("FileArchive lists the following files:")
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
        s = str(type(self)) + "\n"
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

    def read_table_file(self, file, namekey='FILE'):
        """Interprets text in a file input."""

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

    @staticmethod
    def gather_table_from_list(instrument, files):
        """Gather file header information to fill the table

        Args:
            instrument:
            files:

        Returns:

        """

        # Defaults
        header_cards = ['OBSTYPE', 'OBJECT', 'FILTER', 'EXPTIME', 'nFRAMES', 'DATE']
        dtypes = [str, str, str, float, int, str]
        instrument_config_file = os.path.join(os.path.dirname(__file__), '../config/instruments.cfg')

        # Read config
        configs = config.read(instrument_config_file)
        instrument = configs['INSTRUMENTS'][instrument]
        instrument_header_cards = configs[instrument]

        # Double check whether all aliases are defined
        for card in header_cards:
            try:
                instrument_header_cards[card]
            except KeyError:
                logger.info(
                    f"Dropping header card {card} from setup identification, as there is no description in the config "
                    f"file.\nCheck out {instrument_config_file} for details.")
                header_cards.remove(card)

        # Initialize output file information table
        table = Table(names=['FILE']+header_cards, dtype=[str]+dtypes)

        # Read data from files
        for file in files:
            logger.info(f"Retrieving header information from file {file}")
            hdr = fits.getheader(file)
            new_row = [os.path.basename(file)]
            for card in header_cards:
                try:
                    new_row.append(hdr[instrument_header_cards[card]])
                except KeyError:
                    logger.info(f"Skipping file {os.path.basename(file)} due to at least one missing header card "
                                f"({instrument_header_cards[card]}).")
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

    # def update_filenames(self, filenames):
    #     """Updates the file names in the files table.
    #
    #     Args:
    #         filenames (dict):
    #             Dictionary mapping the current file names onto new ones.
    #     """
    #
    #     if not isinstance(filenames, dict):
    #         raise SpecklepyTypeError('FileArchive.update_filenames', argname='filenames', argtype=type(filenames), expected='dict')
    #
    #     logger.info('Updating file names:')
    #     for outdated in filenames.keys():
    #         index = np.where(self.table['FILE'] == outdated)
    #         self.table['FILE'][index] = filenames[outdated]
    #         logger.info(f"{outdated} > {filenames[outdated]}")

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
