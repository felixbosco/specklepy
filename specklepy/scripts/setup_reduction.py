#!/usr/bin/env python

import argparse
import os
import sys
import warnings
import glob
from astropy.io import fits
from astropy.table import Table

try:
    from specklepy.logging import logging
    from specklepy.io.parameterset import ParameterSet
    from specklepy.io.filemanager import FileManager
    from specklepy.reduction.config import get_instrument_config
    from specklepy.reduction import sky
except ModuleNotFoundError:
    # Prepare import from current path
    PATH = os.getcwd()
    warnings.warn("Importing from path {}. Apparently the package is not installed properly on your machine!".format(PATH), ImportWarning)
    sys.path.insert(0, PATH)

    # Repeat import
    from specklepy.logging import logging
    from specklepy.io.parameterset import ParameterSet
    from specklepy.io.filemanager import FileManager
    from specklepy.reduction.config import get_instrument_config
    from specklepy.reduction import sky



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script reduces the data, following the parameters specified in the paramater fils.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--instrument', type=str, help='Name of the instrument.')
    parser.add_argument('-p', '--path', type=str, help='Path to the files.')
    parser.add_argument('-d', '--debug', type=bool, default=False, help='Set to True to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)

    # Default values
    header_cards = ['OBJECT', 'OBSTYPE', 'EXPTIME']
    instrument_config_file = 'specklepy/config/instruments.cfg'

    # Verification of args
    instr = get_instrument_config(args.instrument, instrument_config_file)
    keywords = [instr[card] for card in header_cards]
    if not os.path.isdir(os.path.dirname(args.path)):
        raise RuntimeError("Path not found: {}".format(args.path))

    # Find files
    files = glob.glob(args.path)

    # Prepare table creation with dictionary
    datadict = {'FILE': []}
    for keyword in keywords:
        datadict[keyword] = []

    for file in files:
        hdr = fits.getheader(file)
        try:
            for keyword in keywords:
                datadict[keyword].append(hdr[keyword])
            datadict['FILE'].append(file)
        except KeyError:
            logging.info("Skipping {} due to missing header key ({})".format(file, keyword))

    # Create table from dict and save
    table = Table([datadict[keyword] for keyword in datadict.keys()], names=datadict.keys())
    table.sort('OBSTYPE')
    table.write('files.tab', format='ascii', overwrite=True)







if __name__ == '__main__':
    main()
