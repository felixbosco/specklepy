#!/usr/bin/env python

import argparse
import os
import sys
import warnings

from specklepy.core.holography import holography
from specklepy.io.parameterset import ParameterSet
from specklepy.logging import logging



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script creates a holographic reconstruction of the input files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('parameter_file', type=str, help='Path to the parameter file.')
    parser.add_argument('-m', '--mode', type=str, default='same', help="Reconstruction mode, can be 'same', 'full' or 'valid'.")
    parser.add_argument('-d', '--debug', action='store_true', help='Set to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)

    # Default values
    dir = os.path.dirname(__file__)
    defaults_file = os.path.join(dir, '../config/holography.cfg')
    essential_attributes = {'paths': ['inDir', 'tmpDir', 'outFile', 'alignmentReferenceFile', 'refSourceFile'],
                            'starfinder': ['starfinderFwhm', 'noiseThreshold'],
                            'psfextraction': ['mode', 'psfRadius', 'noiseThreshold', 'noiseReferenceMargin'],
                            'apodization': ['apodizationWidth', 'apodizationType'],
                            'uncertainties': ['varianceExtensionName']}
    make_dirs = ['tmpDir']

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    else:
        params = ParameterSet(parameter_file=args.parameter_file,
                        defaults_file=defaults_file,
                        essential_attributes=essential_attributes,
                        make_dirs=make_dirs)

    # Execute reconstruction
    holography(params, mode=args.mode, debug=args.debug)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logging.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
