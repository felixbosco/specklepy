#!/usr/bin/env python

import argparse
import os
import sys

import numpy as np
from astropy.io import fits
from matplotlib.colors import LogNorm

from specklepy.logging import logger
from specklepy.core.aperture import Aperture
from specklepy.plotting.utils import imshow, plot_simple



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script analyses the statistics of an aperture across the time axis of a cube/ set of cubes for stars in the input file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type=str, default=None, help='File name to analyse.')
    parser.add_argument('-F', '--Fourier_file', type=str, default=None, help='The file may be a tmp file, such that the Fourier transformations do not have to be repeated.')
    parser.add_argument('-i', '--index', nargs='+', type=int, help='Center index of the aperture to analyse. Provide this as "--index 123 456", without other symbols.')
    parser.add_argument('-r', '--radius', type=int, default=10, help='Radius of the aperture to analyse in pix.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the file to write the results to. Default is to just repace .fits by .dat.')
    parser.add_argument('-p', '--pixel_scale', type=float, default=None, help='Pixel scale in arcsec for computing the spatial frequencies.')
    # parser.add_argument('-v', '--visual', action='store_const', const=True, default=False, help='Show the plots?')
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
        if args.Fourier_file is not None:
            start_with_Fourier_file = True
            logger.info("Starting from Fourier file {}.".format(args.Fourier_file))
            args.file = args.Fourier_file
        else:
            raise IOError("At least one of --file or --Fourier_file have to be given!")
    else:
        start_with_Fourier_file = False
        # Apply default
        if args.Fourier_file is None:
            args.Fourier_file = args.file.replace(".fits", "_Fourier.fits")
    # Apply default
    if args.outfile is None:
        args.outfile = args.file.replace(".fits", ".dat")


    # Starting with the pre-analysis of the file
    if not start_with_Fourier_file:

        # Test of the heavy spot of the aperture
        aperture = Aperture(*args.index, args.radius, file_name=args.file, crop=False)
        max = aperture.get_aperture_peak()
        if args.debug:
            imshow(aperture.data, title="Aperture in the integrated cube")
        if max == args.index:
            logger.info("Index {} is identical with the maximum of the integrated cube image {}.".format(args.index, max))
        else:
            logger.info("Index {} is not identical with the maximum of the integrated cube image {}.".format(args.index, max))
            answer = input("Shall the index guess be replaced by the coordinates of the local maximum? (yes,no)")
            if answer.lower() == "yes":
                logger.info("Replacing index {} by the maximum of the integrated cube image {}...".format(args.index, max))
                args.index = max
                aperture = Aperture(*args.index, args.radius, file_name=integrated_cube, mask=None, crop=True)
                if answer and (answer.lower() == "yes" or answer == ''):
                    if args.debug:
                        imshow(aperture.data, title="Updated aperture")
            else:
                logger.info("Continuing with the central index {}...".format(args.index))

        # Remove margins from aperture
        aperture.remove_margins()
        logger.info("The aperture has shape {}.".format(aperture.data.shape))
        aperture.powerspec_to_file(args.file, args.Fourier_file)
        del aperture


    # Show interim results
    Fourier_cube = fits.getdata(args.Fourier_file)
    Fourier_mean = np.mean(Fourier_cube, axis=0)
    Fourier_var = np.var(Fourier_cube, axis=0)  # Compute variance to average linearly in the following
    if args.debug:
        imshow(Fourier_mean, title="Mean of the Fourier transformed cube along the time axis", norm=LogNorm())
        imshow(np.sqrt(Fourier_var), title="Standard deviation in the Fourier transformed cube along the time axis", norm=LogNorm())

    # Compute Fourier radius map and mask
    center = (Fourier_mean.shape[0] / 2, Fourier_mean.shape[1] / 2)
    xx, yy = np.mgrid[:Fourier_mean.shape[0], :Fourier_mean.shape[1]]
    Fourier_radius = np.sqrt(np.square(xx - center[0]) + np.square(yy - center[1]))
    Fourier_radius = np.ma.masked_greater(Fourier_radius, center[0])
    Fourier_mean = np.ma.masked_array(Fourier_mean, mask=Fourier_radius.mask)
    Fourier_var = np.ma.masked_array(Fourier_var, mask=Fourier_radius.mask)

    # Average azimuthally
    logger.info("Averaging the Fourier plane azimuthally")
    Fourier_radius = Fourier_radius.reshape((-1))
    Fourier_mean = Fourier_mean.reshape((-1))
    Fourier_var = Fourier_var.reshape((-1))

    xdata = np.unique(Fourier_radius)
    ydata = np.zeros(xdata.shape)
    edata = np.zeros(xdata.shape)
    for index, value in enumerate(xdata):
        ydata[index] = np.mean(Fourier_mean[np.where(Fourier_radius == value)])
        edata[index] = np.mean(Fourier_var[np.where(Fourier_radius == value)])
    # Turn variance into standard deviation
    edata = np.sqrt(edata)

    if args.pixel_scale is not None:
        logger.warn("Handling the pixel scale is not implemented yet!")

    if args.debug:
        plot_simple(xdata, ydata,
                    title="{}\nCenter={} Radius={}".format(args.Fourier_file, args.index, args.radius),
                    xlabel="Fourier radius")

    # Save power spectra to outfile
    caption = "Fourier_radius mean std"
    data = np.concatenate(([xdata], [ydata], [edata]), axis=0).transpose()
    np.savetxt(args.outfile, data, header=caption)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
