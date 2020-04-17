from configparser import ConfigParser
import glob
import os

from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning

from specklepy.logging import logger


def setup(args):
    header_cards = ['OBSTYPE', 'OBJECT', 'FILTER', 'EXPTIME', 'nFRAMES', 'DATE']
    instrument_config_file = os.path.join(os.path.dirname(__file__), '../config/instruments.cfg')

    # Verification of args
    if not os.path.isdir(os.path.dirname(args.files)):
        raise RuntimeError("Path not found: {}".format(args.files))

    # Read config
    config = ConfigParser()
    config.read(instrument_config_file)
    instrument = config['INSTRUMENTS'][args.instrument]
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
    if '*' in args.files:
        files = glob.glob(args.files)
    else:
        files = glob.glob(args.files + '*fits')
    logger.info("Found {} file(s)".format(len(files)))

    # Prepare dictionary for collecting table data
    table_data = {'FILE': []}
    for card in header_cards:
        table_data[card] = []

    # Read data from files
    for file in files:
        logger.info(f"Retrieving header information from file {file}")
        try:
            hdr = fits.getheader(file)
        except (AstropyWarning, AstropyUserWarning):
            print("Caught")
        table_data['FILE'].append(os.path.basename(file))
        for card in header_cards:
            try:
                table_data[card].append(hdr[instrument_header_cards[card]])
            except KeyError:
                logger.info(
                    f"Skipping file {os.path.basename(file)} due to missing header card ({instrument_header_cards[card]}).")
                table_data[card].append("_" * 3)

    # Create table from dict and save
    table = Table([table_data[keyword] for keyword in table_data.keys()], names=table_data.keys())
    table.sort('FILE')
    table.sort('OBSTYPE')
    table.sort(args.sortby)
    logger.info("Writing data to {}".format(args.outfile))
    table.write(args.outfile, format='ascii.fixed_width', overwrite=True)

    # Write dummy parameter file for the reduction
    logger.info("Creating default reduction INI file {}".format(args.parfile))
    par_file_content = "[PATHS]" \
                       f"\nfilePath = {args.files}" \
                       f"\nfileList = {args.outfile}" \
                       "\ntmpDir = tmp/" \
                       "\n\n[FLAT]" \
                       "\nskipFlat = False" \
                       "\nmasterFlatFile = MasterFlat.fits" \
                       "\nflatCorrectionPrefix = f_" \
                       "\n\n[SKY]" \
                       "\nskipSky = False" \
                       "\nskySubtractionPrefix = s"
    with open(args.parfile, 'w+') as parfile:
        parfile.write(par_file_content)