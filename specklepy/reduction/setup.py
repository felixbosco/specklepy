import glob
from IPython import embed
import numpy as np
import os

from astropy.io import fits
from astropy.table import Column, Table

from specklepy.io import config
from specklepy.logging import logger


def gather_header_information(path, instrument, par_file=None, list_file=None, sort_by=None):
    """Sets up the data reduction parameter file and file list.

    Args:
        path (str):
            Path to the files.
        instrument (str):
            Name of the instrument that took the data. This must be covered by config/instruments.cfg.
        par_file (str, optional):
            Name of the output default parameter file for the reduction.
        list_file (str):
            Name of the output file that contains all the file names and header information.
        sort_by (str, optional):
            Header card that is used for the sorting of files.
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
                f"Dropping header card {card} from setup identification, as there is no description in the config file."
                f"\nCheck out {instrument_config_file} for details.")
            header_cards.remove(card)

    # Apply fall back values
    if path is None:
        path = '.'
    if list_file is None:
        list_file = 'files.tab'
    if par_file is None:
        par_file = 'reduction.yaml'

    # Find files
    if '*' in path:
        files = glob.glob(path)
    else:
        files = glob.glob(os.path.join(path, '*fits'))
    if len(files):
        logger.info(f"Found {len(files)} file(s)")
        files.sort()
    else:
        logger.error(f"Found no files in {path}!")
        raise RuntimeError(f"Found no files in {path}!")

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

    # Sort table entries by default properties and user request
    table.sort('FILE')
    table.sort('OBSTYPE')
    if sort_by:
        table.sort(sort_by)

    # Identify instrument setups
    setups = identify_instrument_setups(table)
    table.add_column(setups)

    # Save table
    logger.info(f"Writing header information to {list_file}")
    table.write(list_file, format='ascii.fixed_width', overwrite=True)

    # Write dummy parameter file for the reduction
    _, ext = os.path.splitext(par_file)
    if 'yaml' in ext:
        logger.info(f"Creating default reduction YAML parameter file {par_file}")
        par_file_content = f"PATHS:\n  filePath: {path}\n  fileList: {list_file}\n  outDir: Science/\n  tmpDir: tmp/" \
                           f"\n\nFLAT:\n  masterFlatFile: MasterFlat.fits" \
                           f"\n\nSKY:\n  method: scalar"
    else:
        logger.info(f"Creating default reduction INI parameter file {par_file}")
        par_file_content = f"[PATHS]\nfilePath = {path}\nfileList = {list_file}\noutDir = Science/\ntmpDir = tmp/" \
                           f"\n\n[FLAT]\nmasterFlatFile = MasterFlat.fits" \
                           f"\n\n[SKY]\nmethod = scalar"
    with open(par_file, 'w+') as par_file:
        par_file.write(par_file_content)


def identify_instrument_setups(table, keys=None):

    # Apply fall back values
    columns = ['FILTER', 'EXPTIME']
    if keys:
        columns = columns + keys

    settings = {}
    n_setups = 1
    for column in columns:
        settings[column] = np.unique(table[column].data)

        # if len(settings[column]) == 1:
        #     columns.pop(column)

        n_setups *= len(settings[column])

    if n_setups == 1:
        return Column(data=np.ones(len(table), dtype=int), name='SETUP', dtype=int)
    else:
        raise NotImplementedError(f"Creating a SETUP table column is not implemented for multiple different setups!")
