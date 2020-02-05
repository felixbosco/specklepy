#!/usr/bin/env python

import argparse
import os
import glob
from configparser import ConfigParser
from astropy.io import fits
from astropy.table import Table

from specklepy.logging import logging



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script reduces the data, following the parameters specified in the paramater fils.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--instrument', type=str, help='Name of the instrument.')
    parser.add_argument('-p', '--path', type=str, help='Path to the files.')
    parser.add_argument('-o', '--outfile', type=str, default='files.tab', help="Name of the output file containing the file overview. Default is 'files.tab'.")
    parser.add_argument('-s', '--sortby', type=str, default='OBSTYPE', help="Header card to sort the output table by. Default is 'OBSTYPE'.")
    parser.add_argument('-d', '--debug', type=bool, default=False, help='Set to True to inspect intermediate results. Default is False.')

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

    # Find files
    files = glob.glob(args.path)
    logging.info("Found {} file(s)".format(len(files)))

    # Prepare dictionary for collecting table data
    table_data = {'FILE': []}
    for card in header_cards:
        table_data[card] = []

    # Read data from files
    for file in files:
        hdr = fits.getheader(file)
        table_data['FILE'].append(os.path.basename(file))
        for card in header_cards:
            try:
                table_data[card].append(hdr[instrument_header_cards[card]])
            except KeyError:
                logging.info("Skipping file {} due to missing header card ({}).".format(os.path.basename(file), instrument_header_cards[card]))
                table_data[card].append("_" * 3)

    # Create table from dict and save
    table = Table([table_data[keyword] for keyword in table_data.keys()], names=table_data.keys())
    table.sort('FILE')
    table.sort('OBSTYPE')
    table.sort(args.sortby)
    logging.info("Writing data to {}".format(args.outfile))
    table.write(args.outfile, format='ascii.fixed_width', overwrite=True)



if __name__ == '__main__':
    main()
