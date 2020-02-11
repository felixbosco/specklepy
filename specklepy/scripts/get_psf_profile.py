#!/usr/bin/env python

import argparse
import os
import sys
import warnings
import numpy as np

from specklepy.core.aperture import Aperture
from specklepy.utils.plot import imshow, encircled_energy_plot



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script extracts the encircled energy as a function of readius within an aperture.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='File name to analyse.')
    parser.add_argument('-i', '--index', nargs='+', type=int, help='Center index of the aperture to analyse. Provide this as "--index 123 456", without other symbols.')
    parser.add_argument('-r', '--radius', type=int, default=10, help='Radius of the aperture to analyse in pix.')
    parser.add_argument('-n', '--normalize', type=str, default=None, help='Normalize the flux values, to either "peak", "aperture" or set to None for not normalizing. Default is None.')
    parser.add_argument('-o', '--outdir', type=str, default=None, help='Directory to write the results to. Default is to repace .fits by .dat and add a "energy" prefix.')
    parser.add_argument('-m', '--maximize', type=bool, default=True, help='Set to True to show every debug plot on full screen. Default is True.')
    parser.add_argument('-d', '--debug', type=bool, default=False, help='Set to True to inspect results.')

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

    outfile = os.path.basename(args.file)
    outfile = "psf_" + outfile.replace(".fits", ".dat")
    outfile = os.path.join(args.outdir, outfile)


    # SInitialize the aperture
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
        encircled_energy_plot(outfile, maximize=args.maximize)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logging.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
