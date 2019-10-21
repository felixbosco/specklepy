import numpy as np
import argparse
from astropy.io import fits
# from matplotlib.colors import LogNorm

try:
    from holopy.logging import logging
    from holopy.core.aperture import Aperture
    from holopy.utils.plot import imshow
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
    from holopy.utils.plot import imshow


def parser(options=None):

    parser = argparse.ArgumentParser(description='This script analyses the statistics of an aperture across the time axis of a cube/ set of cubes for stars in the input file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='Generic file name to analyse.')
    parser.add_argument('-i', '--index', nargs='+', type=int, help='Center index of the aperture to analyse. Provide this as "--index 123 456", without other symbols.')
    parser.add_argument('-r', '--radius', type=int, default=10, help='Radius of the aperture to analyse in pix.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(options=None):

    args = parser(options=options)

    # Prepare input
    args.index = tuple(args.index)
    args.tmpfile = args.file.replace(".fits", "_tmp.fits")
    args.outfile = args.file.replace(".fits", ".dat")
    # print(args)

    # Load cube
    cube = fits.getdata(args.file)

    # Test of the heavy spot of the aperture
    moment0 = np.sum(cube, axis=0)
    aperture = Aperture(*args.index, args.radius, data=moment0, subset_only=False)
    max = np.unravel_index(np.argmax(aperture(), axis=None), aperture().shape)
    imshow(aperture())
    if max == args.index:
        logging.info("Index {} is identical with the maximum of the integrated cube image {}.".format(args.index, max))
    else:
        logging.info("Index {} is not identical with the maximum of the integrated cube image {}.".format(args.index, max))
        answer = input("Shall the index guess be replaced by the coordinates of the local maximum? (yes,no)")
        if answer.lower() == "yes":
            logging.info("Replacing index {} by the maximum of the integrated cube image {}...".format(args.index, max))
            args.index = max
        else:
            logging.info("Continuing with the central index {}...".format(args.index))
    del aperture

    # Analysis
    aperture = Aperture(*args.index, args.radius, data=moment0, mask=None, subset_only=True)
    imshow(aperture.data)




if __name__ == '__main__':
    main()
