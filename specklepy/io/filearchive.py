from datetime import datetime
import glob
import numpy as np
import os
import string
import sys

from astropy.io.registry import IORegistryError
from astropy.table import Column, Table

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.fits import get_header
from specklepy.io.filestream import FileStream
from specklepy.logging import logger
from specklepy.reduction.sequence import Sequence
from specklepy.utils.time import default_time_stamp


class FileArchive(object):

    """This class is handling the input files and is iterable.

    Attributes:
        file_list (str):
            Path to list of files or generic file path.
        ...
    """

    def __init__(self, file_list, file_path=None, out_dir=None, out_prefix=None, cards=None, dtypes=None,
                 table_format=None):
        """Create a FileArchive instance.

        Long description...

        Args:
            file_list (str, list):
                Path to list of files or generic file path. Can also be provided as list type.
            file_path (str, optional):
                Path to the raw/ input data.
            out_dir (str, optional):
                Path to the product/ output data.
            out_prefix (str, optional):
                Prefix of the product/ output data.
            cards (list, optional):
                Header keywords, for which the values shall be extracted from the FITS headers.
            dtypes (list, optional):
                List of data types, to which the header entries of the `cards` are casted.
            table_format (str, optional):
                Format of the ASCII file, from which file names are read. Used only if `file_list` is parsed as the name
                of an ASCII table file.
        """

        # Store in and out paths
        self.file_path = file_path
        self.out_dir = out_dir if out_dir is not None else './'
        self.out_prefix = out_prefix if out_prefix is not None else ''

        # Interpret the str-type file list input
        if isinstance(file_list, str):
            file_list = self.find_files(file_list, path=self.file_path)

        # Assert type of file list
        if not isinstance(file_list, list):
            raise SpecklepyTypeError("FileArchive", 'file_list', type(file_list), 'list')

        # Read file names from list file
        if len(file_list) == 1 and not self.is_fits_file(file_list[0]):
            logger.info(f"Input file is not FITS type. FileArchive assumes that input file {file_list[0]!r} contains "
                        f"file names.")
            file_list = self.read_table_file(file_list[0], format=table_format)

        # Gather table of file names and request further information, from the FITS header
        self.gather_table_from_list(files=file_list, cards=cards, dtypes=dtypes)

        # Log identified input files
        logger.debug("FileArchive lists the following files:")
        logger.debug(str(self.table))

        # Initialize the index for iteration
        self.index = 0

        # Reduce paths
        self.file_path = os.path.normpath(self.file_path)
        self.out_dir = os.path.normpath(self.out_dir)

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

    @staticmethod
    def find_files(generic, path=None):

        # Join path into file name
        if path is not None:
            generic = os.path.join(path, generic)

        files = glob.glob(generic)
        files.sort()
        if len(files) == 0:
            sys.tracebacklimit = 0
            raise FileNotFoundError(f"FileArchive did not find any file matching to {generic!r}.")
        else:
            logger.info(f"FileArchive found {len(files)} file(s) matching to {generic!r}.")

        return files

    def read_table_file(self, file, format='ascii.no_header'):
        """Interprets text in a file input.

        Args:
            file (str):
                Name of the file containing the table.
            format (str, optional):
                Format string of the input table.

        Returns:
            table (astropy.Table):
                Table based on file input.
            common_path (str):
                Common path to the files in the table.
        """

        logger.info(f"Reading file names from input file {file!r}")

        # Read table from ASCII file
        try:
            table = Table.read(file, format=format)
            if format == 'ascii.no_header':
                default_column_name = table.colnames[0]
                table[default_column_name].name = 'FILE'
        except IORegistryError as e:
            sys.tracebacklimit = 0
            raise e

        # Return list of files
        return list(table['FILE'].data)

    def gather_table_from_list(self, files, cards=None, dtypes=None, names=None, sort_by=None):
        """Gather file header information to fill the table

        Args:
            files (list):
                Names of the FITS files to read the header information from.
            cards (list):
                Header cards to read.
            dtypes (list):
                Expected data types of values.
            names (list):
                Names of the output table columns. If provided this overwrites the column names based on `cards`.
            sort_by (str, optional):
                Header card that is used for the sorting of files.

        Returns:
            table (astropy.Table):
                Table with the header information about the header `cards`, for the files in `files`.
        """

        # Initialize output file information table
        if cards is None:
            cards = []
        if dtypes is None:
            dtypes = []
        if names is None:
            table = Table(names=['FILE']+cards, dtype=[str]+dtypes)
        else:
            table = Table(names=['FILE']+names, dtype=[str]+dtypes)

        # Extract common path
        self.file_path = self.common_path(files)

        # Read data from files
        for path in files:
            file = path.replace(self.file_path, '').lstrip('/')

            logger.info(f"Retrieving header information from file {file!r}")
            hdr = get_header(file, path=self.file_path)
            new_row = [file]
            for card in cards:
                # Card actually contains two or more cards
                if ',' in card:
                    _cards = card.replace(',', '').split()
                    value = hdr[_cards[0]]
                    for _card in _cards[1:]:
                        value += ' ' + hdr[_card]
                    new_row.append(value)
                # Only one card
                else:
                    if card in hdr:
                        value = hdr[card]
                        if not isinstance(value, str) or card == 'FILE':
                            new_row.append(value)
                        else:
                            new_row.append(value.upper())
                    elif card == 'NAXIS3':
                        logger.warning(f"File {file!r} is missing header card {card!r}. Number of frames set to 1!")
                        new_row.append(1)
                    else:
                        new_row.append(None)
            if len(new_row) == len(table.columns):
                table.add_row(new_row)

        # Sort table entries by default properties and user request
        try:
            table.sort('DATE')
        except ValueError:
            pass  # 'DATE' is not a valid column in table
        if sort_by:
            table.sort(sort_by)

        self.table = table

    @staticmethod
    def common_path(files):
        if len(files) == 1:
            return os.path.dirname(files[0])
        else:
            return os.path.commonpath(files)

    def write_table(self, file_name):
        """Write the archive's table to a file `file_name`."""
        logger.info(f"Writing file information to {file_name!r}")
        self.table.write(file_name, format='ascii.fixed_width', overwrite=True)


class ReductionFileArchive(FileArchive):

    @property
    def setups(self):
        return np.unique(self.table['SETUP'].data)

    @property
    def objects(self):
        return np.unique(self.table['OBJECT'].data)

    @property
    def source_files(self):
        return self.table['FILE'].data

    @property
    def source_file_paths(self):
        paths = []
        for source_file in self.source_files:
            paths.append(os.path.join(self.file_path, source_file))
        return paths

    @property
    def product_files(self):
        try:
            out =  self.table['PRODUCT'].data
        except KeyError:
            self.add_product_file_column()
            out = self.table['PRODUCT'].data
        return out

    @property
    def product_file_paths(self):
        paths = []
        for product_file in self.product_files:
            if product_file is not None:
                product_file_path = os.path.join(self.out_dir, product_file)
            else:
                product_file_path = None
            paths.append(product_file_path)
        return paths

    def filter(self, filter_dict, namekey='FILE', return_mask=False):
        """Filter the archive's table by the column properties.

        filter_dict = {'name_of_column': [desired_value_1, desired_value_2]}

        Args:
            filter_dict:
                Dictionary that holds the table column names as keys and acceptable values as lists as the
                corresponding values.
            namekey (str, optional):
                Name/ key of the out put column that shall be filtered by `filter_dict`.
            return_mask (bool, optional):
                Return also the full mask.

        Returns:
            filtered (np.array):
                Array_like subset of the column `namekey` that match the criteria based on filter_dict.
        """

        if not isinstance(filter_dict, dict):
            raise SpecklepyTypeError('FileArchive.filter', 'filter_dict', type(filter_dict), 'dict')

        mask = [True] * len(self.table[namekey])
        for key, search_value in filter_dict.items():
            if isinstance(search_value, list):
                submask = [False] * len(self.table[namekey])
                for correct in search_value:
                    submask |= (self.table[key] == correct)
                mask &= submask
            else:
                mask &= (search_value == self.table[key])

        if return_mask:
            return self.table[namekey][mask].data, mask

        return self.table[namekey][mask].data

    def get_dark_setups(self):
        is_dark = self.table['OBSTYPE'] == 'DARK'
        return np.unique(self.table['SETUP'][is_dark].data)

    def add_dark_column(self):
        """Add a column to the table and initialize with dark setups of matching exposure time."""

        # Initialize column
        dark_col = Column(name='DARK', dtype=object, length=len(self.table))

        # Fill the column with the dark setup of matching exposure time
        dark_setups = self.get_dark_setups()
        for setup in dark_setups:
            exp_times = self.filter({'OBSTYPE': 'DARK', 'SETUP': setup}, namekey='EXPTIME')
            exp_time = np.unique(exp_times)[0]
            dark_col[self.table['EXPTIME'] == exp_time] = setup

        # Append new column to the table
        self.table.add_column(dark_col)

    def get_flats(self):
        return self.filter({'OBSTYPE': 'FLAT'})

    def get_flat_filters(self):
        is_flat = self.table['OBSTYPE'] == 'FLAT'
        return np.unique(self.table['FILTER'][is_flat].data)

    def get_science(self):
        return self.filter({'OBSTYPE': 'SCIENCE'})

    def get_skies(self):
        return self.filter({'OBSTYPE': 'SKY'})

    def identify_setups(self, keywords):
        """Identify distinct observational setups in a list of files.

        This function identifies distinct observational setups.

        Args:
            keywords (list of str):
                Keywords that shall be used to distinguish between observational setups.
        """

        # Check input parameters
        if not isinstance(keywords, list):
            raise SpecklepyTypeError('identify_setups', 'keywords', type(keywords), 'list')

        # Identifying setups key-by-key
        logger.info("Identifying distinct observational setups in the file list...")
        self.table.add_column(col=Column(data=[None] * len(self.table), name='SETUP'))

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
                            row['SETUP'] += str(index)

        # Overwrite setup keys by length-1 string
        combinations = np.unique(self.table['SETUP'].data)
        for index, combination in enumerate(combinations):
            row_indexes = np.where(self.table['SETUP'].data == combination)
            self.table['SETUP'][row_indexes] = string.ascii_uppercase[index]

    def identify_sequences(self, source='sky'):
        """Identify observation sequences.

        Args:
            source (str, optional):
                Observation type of the images the shall be used to measure the sky background from. Options are 'sky'
                (default) and 'science'.

        Returns:
            sequences (list of Sequence):
                List of observing sequences.
        """

        # Type check
        if isinstance(source, str):
            if source not in ['sky', 'science']:
                raise SpecklepyValueError('identify sequences', argname='source', argvalue=source,
                                          expected="'sky' or 'science'")
        else:
            raise SpecklepyTypeError('identify sequences', argname='source', argtype=type(source), expected='str')

        # Identify the observing sequences
        sequences = []
        for setup in self.setups:
            for object in self.objects:
                # Query names and time stamps of science and sky files
                sky_files = self.filter({'OBSTYPE': source.upper(), 'OBJECT': object, 'SETUP': setup})
                sky_time_stamps = self.filter({'OBSTYPE': source.upper(), 'OBJECT': object, 'SETUP': setup},
                                              namekey='DATE')
                science_files = self.filter({'OBSTYPE': 'SCIENCE', 'OBJECT': object, 'SETUP': setup})
                science_time_stamps = self.filter({'OBSTYPE': 'SCIENCE', 'OBJECT': object, 'SETUP': setup},
                                                  namekey='DATE')

                # Test the number of source files
                if len(sky_files) == 0:
                    logger.warning(f"Did not find any sky observations for object {object} in setup {setup}. No sky "
                                   f"subtraction will be applied!")
                else:
                    # Store the information in a new sequence
                    sequences.append(Sequence(sky_files=sky_files, science_files=science_files, file_path=self.file_path,
                                              sky_time_stamps=sky_time_stamps, science_time_stamps=science_time_stamps,
                                              source=source, object=object, setup=setup))
        return sequences

    def add_product_file_column(self, prefix=None):
        """Add a column of product file names to the table.

        Arguments:
            prefix (str, optional):
                .
        """

        # Store update prefix
        if prefix:
            self.out_prefix = prefix

        # Initialize column
        product_file_column = Column(name='PRODUCT', dtype=object, length=len(self.table))

        # Define list of obs-types, considered for product files
        obs_types = ['SCIENCE', 'SKY', 'FLAT']

        # Iterate through table and fill column
        for r, row in enumerate(self.table):
            if row['OBSTYPE'] in obs_types:
                dest = self.out_prefix + os.path.basename(row['FILE'])
                product_file_column[r] = dest
            else:
                logger.debug(f"Ignoring file {row['FILE']!r} for its OBSTYPE {row['OBSTYPE']!r} not being in the list "
                             f"of considered types: {obs_types}")

        # Fill empty fields
        is_empty = product_file_column == 0
        product_file_column[is_empty] = None

        # Store the new column of product file names to the table
        self.table.add_column(product_file_column)

    def initialize_product_file(self, index, prefix=None):
        """Copy the science data cubes into the stored out directory.

        Args:
            index (int):
                Index of the file to be initialized in the list of `self.source_files`.
            prefix (str, optional):
                File prefix for output files.

        Returns:
            product_files (list):
                List of paths of the data reduction products.
        """

        # Copy the science data cubes into outdir (with an additional file prefix)
        src = self.source_file_paths[index]
        dest = self.product_file_paths[index]
        logger.info(f"Initializing data product file {dest!r}")
        os.system(f"cp {src} {dest}")

        # Add information to the header of the new FITS file
        product_file = FileStream(file_name=dest)
        product_file.set_header('PIPELINE', 'SPECKLEPY')
        product_file.set_header('REDUCED', default_time_stamp())

        return dest
