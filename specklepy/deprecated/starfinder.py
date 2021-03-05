#!/usr/bin/env python

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture

from specklepy.core.sourceextraction import extract_sources
from specklepy.io.parameterset import ParameterSet
from specklepy.logging import logger



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script searches for stars in the input file and creates a list file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type=str, default=None, help='Fits file to search for stars.')
    parser.add_argument('-p', '--parameter_file', type=str, help='Path to the parameter file.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the output file.')
    parser.add_argument('-d', '--debug', action='store_true', help='Set to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)

    # Default values
    background_subtraction = True
    defaults_file = os.path.join(os.path.dirname(__file__), "../config/holography.cfg")
    essential_attributes = {'paths': ['allStarsFile'],
                            # 'noise': ['noiseBoxX', 'noiseBoxY', 'noiseBoxHalfWidth'],
                            'starfinder': ['noiseThreshold', 'starfinderFwhm']}

    if args.file is None:
        logger.error("No file or file list was provided! Use --help for instructions.")
        raise RuntimeError("No file or file list was provided! Use --help for instructions.")

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    else:
        params = ParameterSet(parameter_file=args.parameter_file, defaults_file=defaults_file, essential_attributes=essential_attributes)

    if args.outfile is None:
        args.outfile = params.paths.allStarsFile

    # finder = SourceExtraction()
    # finder.
    extract_sources(image=args.file, fwhm=params.starfinder.starfinderFwhm, noise_threshold=params.starfinder.noiseThreshold,
                    background_subtraction=background_subtraction, write_to=args.outfile)
    # finder.writeto(args.outfile)

    # Plot results
    if args.debug:
        positions = np.transpose((finder.sources['xcentroid'], finder.sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.)
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(image, cmap='Greys', origin='lower', norm=norm)
        plt.colorbar(pad=0.0)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.show()
        plt.close()



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
