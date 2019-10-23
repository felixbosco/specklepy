import numpy as np
import argparse
from astropy.io import fits
from datetime import datetime
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

try:
    from holopy.logging import logging
    from holopy.core.aperture import Aperture
    from holopy.utils.plot import imshow
    from holopy.utils import transferfunctions as tf
except ModuleNotFoundError:
    # Prepare import with hardcoded path
    import warnings
    PATH = '/home/bosco/Documents/phd/sowat/pipeline/github_holopy'
    warnings.warn("Importing holopy from hardcoded path {}. Apparently holopy is not installed properly on your machine!".format(PATH), ImportWarning)
    import sys
    sys.path.insert(0, PATH)

    # Repeat import
    from holopy.logging import logging
    from holopy.core.aperture import Aperture
    from holopy.utils.plot import imshow, plot_simple
    from holopy.utils import transferfunctions as tf


def parser(options=None):

    parser = argparse.ArgumentParser(description='This script analyses the statistics of an aperture across the time axis of a cube/ set of cubes for stars in the input file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='File name to analyse.')
    parser.add_argument('-F', '--Fourier_file', type=str, default=None, help='The file may be a tmp file, such that the Fourier transformations do not have to be repeated.')
    parser.add_argument('-i', '--index', nargs='+', type=int, help='Center index of the aperture to analyse. Provide this as "--index 123 456", without other symbols.')
    parser.add_argument('-r', '--radius', type=int, default=10, help='Radius of the aperture to analyse in pix.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the file to write the results to. Default is to just repace .fits by .dat.')
    parser.add_argument('-p', '--pixel_scale', type=float, default=None, help='Pixel scale in arcsec for computing the spatial frequencies.')
    parser.add_argument('-v', '--visual', type=bool, default=False, help='Show the plots?')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(options=None):

    args = parser(options=options)


    # Interprete input
    args.index = tuple(args.index)
    if args.visual:
        raise ValueError('')
    if args.file is None:
        if args.Fourier_file is not None:
            start_with_Fourier_file = True
            logging.info("Starting from Fourier file {}.".format(args.Fourier_file))
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
        cube = fits.getdata(args.file)

        # Test of the heavy spot of the aperture
        moment0 = np.sum(cube, axis=0)
        aperture = Aperture(*args.index, args.radius, data=moment0, subset_only=False)
        max = np.unravel_index(np.argmax(aperture.data, axis=None), aperture.data.shape)
        if args.visual:
            imshow(aperture.data, title="Aperture in the integrated cube")
        if max == args.index:
            logging.info("Index {} is identical with the maximum of the integrated cube image {}.".format(args.index, max))
        else:
            logging.info("Index {} is not identical with the maximum of the integrated cube image {}.".format(args.index, max))
            answer = input("Shall the index guess be replaced by the coordinates of the local maximum? (yes,no)")
            if answer.lower() == "yes":
                logging.info("Replacing index {} by the maximum of the integrated cube image {}...".format(args.index, max))
                args.index = max
                aperture = Aperture(*args.index, args.radius, data=moment0, mask=None, subset_only=True)
                if answer and (answer.lower() == "yes" or answer == ''):
                    if args.visual:
                        imshow(aperture.data, title="Updated aperture")
            else:
                logging.info("Continuing with the central index {}...".format(args.index))
        del aperture

        # Redefine aperture with the updated index and remove margins
        aperture = Aperture(*args.index, args.radius, data=cube, mask=None, subset_only=True)
        del cube
        logging.info("The aperture has shape {}.".format(aperture.data.shape))

        # Initialize Fourier File
        logging.info("Initializing Fourier file {}".format(args.Fourier_file))
        header = fits.getheader(args.file)
        header.set('HIERARCH HOLOPY TYPE', 'Fourier transform of an aperture')
        header.set('HIERARCH HOLOPY ORIGIN', args.file)
        header.set('HIERARCH HOLOPY APERTURE INDEX', str(args.index))
        header.set('HIERARCH HOLOPY APERTURE RADIUS', args.radius)
        header.set('UPDATED', str(datetime.now()))
        data = np.zeros(aperture.data.shape)
        fits.writeto(args.Fourier_file, data=data, header=header, overwrite=True)
        logging.info("Initialized {}".format(args.Fourier_file))

        # Fourier transform analysis
        with fits.open(args.Fourier_file, mode='update') as hdulist:
            for index, frame in enumerate(aperture.data):
                print("\rFourier transforming frame {}/{}".format(index+1, aperture.data.shape[0]), end='')
                hdulist[0].data[index] = tf.powerspec(frame)
                hdulist.flush()
            print()
        logging.info("Computed the Fourier transform of every frame and saved them to {}".format(args.Fourier_file))
        del aperture


    # Show interim results
    Fourier_cube = fits.getdata(args.Fourier_file)
    Fourier_mean = np.mean(Fourier_cube, axis=0)
    Fourier_var = np.var(Fourier_cube, axis=0)  # Compute variance to average linearly in the following
    if args.visual:
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
    logging.info("Averaging the Fourier plane azimuthally")
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
        logging.warn("Handling the pixel scale is not implemented yet!")

    if args.visual:
        plot_simple(xdata, ydata,
                    title="{}\nCenter={} Radius={}".format(args.Fourier_file, args.index, args.radius),
                    xlabel="Fourier radius")

    # Save power spectra to outfile
    caption = "Fourier_radius mean std"
    data = np.concatenate(([xdata], [ydata], [edata]), axis=0).transpose()
    np.savetxt(args.outfile, data, header=caption)









if __name__ == '__main__':
    main()
