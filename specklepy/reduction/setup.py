from configparser import ConfigParser
import glob
import os

from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning

from specklepy.logging import logger


def setup(path, instrument, parfile, filelist, sortby=None):
    """Sets up the data reduction parameter file and file list.

    Args:
        path (str):
            Path to the files.
        instrument (str):
            Name of the instrument that took the data. This must be covered by config/instruments.cfg.
        parfile (str):
            Name of the parameter file.
        filelist (str):
            Name of the file that contains all the files.
        sortby (str, optional):
            Header card that is used for the sorting of files.
    """

    # Defaults
    header_cards = ['OBSTYPE', 'OBJECT', 'FILTER', 'EXPTIME', 'nFRAMES', 'DATE']
    instrument_config_file = os.path.join(os.path.dirname(__file__), '../config/instruments.cfg')

    # Verification of args
    if path is None:
        path = '.'
    if not os.path.isdir(path):
        raise RuntimeError(f"Path not found: {path}")

    # Read config
    config = ConfigParser()
    config.read(instrument_config_file)
    instrument = config['INSTRUMENTS'][instrument]
    instrument_header_cards = config[instrument]

    # Double check whether all aliases are defined
    for card in header_cards:
        try:
            instrument_header_cards[card]
        except:
            logger.info(
                f"Dropping header card {card} from setup identification, as there is no description in the config file."
                f"\nCheck out {instrument_config_file} for details.")
            header_cards.remove(card)

    # Find files
    if '*' in path:
        files = glob.glob(path)
    else:
        files = glob.glob(os.path.join(path, '*fits'))
    if len(files):
        logger.info(f"Found {len(files)} file(s)")
    else:
        logger.error(f"Found no files in {path}!")
        raise RuntimeError(f"Found no files in {path}!")

    # Prepare dictionary for collecting table data
    table_data = {'FILE': []}
    for card in header_cards:
        table_data[card] = []

    # Read data from files
    for file in files:
        logger.info(f"Retrieving header information from file {file}")
        # try:
        hdr = fits.getheader(file)
        # except (AstropyWarning, AstropyUserWarning):
        #     print("Caught")
        table_data['FILE'].append(os.path.basename(file))
        for card in header_cards:
            try:
                table_data[card].append(hdr[instrument_header_cards[card]])
            except KeyError:
                logger.info(f"Skipping file {os.path.basename(file)} due to missing header card "
                            f"({instrument_header_cards[card]}).")
                table_data[card].append("/" * 3)

    # Create table from dict and save
    table = Table([table_data[keyword] for keyword in table_data.keys()], names=table_data.keys())
    table.sort('FILE')
    table.sort('OBSTYPE')
    if sortby:
        table.sort(sortby)
    logger.info(f"Writing header information to {filelist}")
    table.write(filelist, format='ascii.fixed_width', overwrite=True)

    # Write dummy parameter file for the reduction
    logger.info(f"Creating default reduction INI file {parfile}")
    par_file_content = "[PATHS]" \
                       f"\nfilePath = {path}" \
                       f"\nfileList = {filelist}" \
                       "\ntmpDir = tmp/" \
                       "\n\n[FLAT]" \
                       "\nmasterFlatFile = MasterFlat.fits" \
                       "\nflatCorrectionPrefix = f_" \
                       "\n\n[SKY]" \
                       "\nskySubtractionPrefix = s"
    with open(parfile, 'w+') as parfile:
        parfile.write(par_file_content)