#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np

from specklepy.core.aperture import Aperture
from specklepy.logging import logger
from specklepy.plotting.utils import imshow, psf_profile_plot



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script extracts the PSF profile as a function of radius within an aperture.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type=str, default=None, help='File name to analyse.')
    parser.add_argument('-i', '--index', nargs='+', type=int, help='Center index of the aperture to analyse. Provide this as "--index 123 456", without other symbols.')
    parser.add_argument('-r', '--radius', type=int, default=10, help='Radius of the aperture to analyse in pix.')
    parser.add_argument('-n', '--normalize', type=str, default=None, help='Normalize the flux values, to either "peak", "aperture" or set to None for not normalizing. Default is None.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Saves the results to this file.')
    parser.add_argument('-m', '--maximize', action='store_true', help='Set to show every debug plot on full screen.')
    parser.add_argument('-d', '--debug', action='store_true', help='Set to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)


    # Interprete input
    args.index = tuple(args.index)

    if args.file is None:
        raise RuntimeError("No file was provided!")

    if args.outfile is None:
        outfile = os.path.basename(args.file)
        outfile = "psf_" + outfile.replace(".fits", ".dat")
        # outfile = os.path.join(args.outdir, outfile)


    # Initialize the aperture
    aperture = Aperture(args.index, args.radius, data=args.file, crop=True)
    # peak = aperture.get_aperture_peak()
    # aperture = Aperture(peak, args.radius, data=args.file, crop=True)
    if args.debug:
        imshow(aperture.get_integrated(), maximize=args.maximize)


    xdata, ydata = aperture.get_psf_profile()

    if args.normalize == 'peak':
        ydata /= ydata[0]
    elif args.normalize == 'aperture':
        ydata /= ydata[-1]
    elif args.normalize is not None:
        raise ValueError("Normalize must be either 'peak', 'aperture, or None!'")


    # Save encircled energy data to outfile
    header = "Radius Flux"
    data = np.concatenate(([xdata], [ydata]), axis=0).transpose()
    np.savetxt(outfile, data, header=header)

    if args.debug:
        psf_profile_plot(outfile, maximize=args.maximize)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
