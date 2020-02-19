#!/usr/bin/env python

import argparse
import os
import sys
import warnings

from specklepy.logging import logging
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script generates synthetic exposures from a parameter file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('parameter_file', type=str, help='Path to the parameter file.')
    parser.add_argument('-d', '--debug', action='store_true', help='Set to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    objects = get_objects(args.parameter_file, debug=args.debug)

    # Generate exposures
    generate_exposure(target=objects['target'],
                        telescope=objects['telescope'],
                        detector=objects['detector'],
                        **objects['kwargs'],
                        debug=args.debug)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logging.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
