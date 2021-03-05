#!/usr/bin/env python

import argparse
import os
import sys
import glob
import warnings
from configparser import ConfigParser
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning

from specklepy.logging import logger



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script sets up the file list and paramater file for data reduction.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('path', type=str, help='Path to the files.')
    parser.add_argument('-i', '--instrument', type=str, help='Name of the instrument.')
    parser.add_argument('-o', '--outfile', type=str, default='specklepy_reduction_files.tab', help="Name of the output file containing the file overview.")
    parser.add_argument('-p', '--parfile', type=str, default='specklepy_reduction.ini', help="Name of the output parameter file.")
    parser.add_argument('-s', '--sortby', type=str, default='OBSTYPE', help="Header card to sort the output table by.")
    parser.add_argument('-d', '--debug', action='store_true', help='Set to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)

    # Default values
    header_cards = ['OBSTYPE', 'OBJECT', 'FILTER', 'EXPTIME', 'nFRAMES', 'DATE']
    instrument_config_file = os.path.join(os.path.dirname(__file__), '../config/instruments.cfg')

    # Verification of args
    if not os.path.isdir(os.path.dirname(args.path)):
        raise RuntimeError("Path not found: {}".format(args.path))

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
            logger.info(f"Dropping header card {card} from setup identification, as there is no description in the config file."
                         f"\nCheck out {instrument_config_file} for details.")
            header_cards.remove(card)

    # Find files
    if '*' in args.path:
        files = glob.glob(args.path)
    else:
        files = glob.glob(args.path + '*fits')
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
                logger.info(f"Skipping file {os.path.basename(file)} due to missing header card ({instrument_header_cards[card]}).")
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
    par_file_content = "[PATHS]"\
                       f"\nfilePath = {args.path}"\
                       f"\nfileList = {args.outfile}"\
                       "\ntmpDir = tmp/"\
                       "\n\n[FLAT]"\
                       "\nskipFlat = False"\
                       "\nmasterFlatFile = MasterFlat.fits"\
                       "\nflatCorrectionPrefix = f_"\
                       "\n\n[SKY]"\
                       "\nskipSky = False"\
                       "\nskySubtractionPrefix = s"
    with open(args.parfile, 'w+') as parfile:
        parfile.write(par_file_content)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
