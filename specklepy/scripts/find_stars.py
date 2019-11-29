#!/usr/bin/env python

import argparse
import os
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture


try:
    from specklepy.logging import logging
    from specklepy.io.parameterset import ParameterSet
    from specklepy.core.aperture import Aperture
    from specklepy.utils.plot import imshow
    from specklepy.algorithms.sourceextraction import SourceExtraction
except ModuleNotFoundError:
    # Prepare import with hardcoded path
    import warnings
    PATH = os.getcwd()
    warnings.warn("Importing specklepy from hardcoded path {}. Apparently specklepy is not installed properly on your machine!".format(PATH), ImportWarning)
    import sys
    sys.path.insert(0, PATH)

    # Repeat import
    from specklepy.logging import logging
    from specklepy.io.parameterset import ParameterSet
    from specklepy.core.aperture import Aperture
    from specklepy.utils.plot import imshow
    from specklepy.algorithms.sourceextraction import SourceExtraction


def parser(options=None):

    parser = argparse.ArgumentParser(description='This script searches for stars in the input file and creates a list file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='Fits file to search for stars.')
    parser.add_argument('-p', '--parameter_file', type=str, help='Path to the parameter file.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the output file.')
    parser.add_argument('-d', '--debug', type=bool, default=False, help='Set to True to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(options=None):

    args = parser(options=options)

    # Default values
    background_subtraction = True
    defaults_file = "specklepy/config/holography.cfg"
    essential_attributes = ['allStarsFile', 'noiseBoxX', 'noiseBoxY', 'noiseBoxHalfWidth', 'noiseThreshold', 'starfinderFwhm']

    if args.file is None:
        logging.error("No file or file list was provided! Use --help for instructions.")
        raise RuntimeError("No file or file list was provided! Use --help for instructions.")

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    else:
        params = ParameterSet(parameter_file=args.parameter_file, defaults_file=defaults_file, essential_attributes=essential_attributes)

    if args.outfile is None:
        args.outfile = params.allStarsFile

    finder = SourceExtraction()
    finder.find_sources(image=args.file, starfinder_fwhm=params.starfinderFwhm, noise_threshold=params.noiseThreshold,
        background_subtraction=background_subtraction)
    finder.writeto(args.outfile)

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
    main()
