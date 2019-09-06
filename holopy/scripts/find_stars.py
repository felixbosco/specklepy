import argparse
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture


try:
    from holopy.logging import logging
    from holopy.io.paramhandler import ParamHandler
    from holopy.core.aperture import Aperture
    from holopy.utils.imshow import imshow
except ModuleNotFoundError:
    # Prepare import with hardcoded path
    import warnings
    PATH = '/home/bosco/Documents/phd/sowat/pipeline/github_holopy'
    warnings.warn("Importing holopy from hardcoded path {}. Apparently holopy is not installed properly on your machine!".format(PATH), ImportWarning)
    import sys
    sys.path.insert(0, PATH)

    # Repeat import
    from holopy.logging import logging
    from holopy.io.paramhandler import ParamHandler
    from holopy.core.aperture import Aperture
    from holopy.utils.imshow import imshow


def parser(options=None):

    parser = argparse.ArgumentParser(description='This script searches for stars in the input file and creates a list file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='Fits file to search for stars.')
    parser.add_argument('-p', '--parameter_file', type=str, help='Path to the parameter file.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the output file.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(options=None):

    args = parser(options=options)

    # Default values
    plot = False
    background_subtraction = True
    defaults_file = "holopy/config/holography_defaults.cfg"
    essential_attributes = ['allStarsFile', 'noiseBoxX', 'noiseBoxY', 'noiseBoxHalfWidth', 'noiseThreshold', 'starfinderFwhm']

    if args.file is None:
        logging.error("No file or file list was provided! Use --help for instructions.")
        raise RuntimeError("No file or file list was provided! Use --help for instructions.")

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    else:
        params = ParamHandler(parameter_file=args.parameter_file, defaults_file=defaults_file, essential_attributes=essential_attributes)

    # Prepare noise statistics
    image = fits.getdata(args.file)
    # noise_box = Aperture(params.noiseBoxX, params.noiseBoxY, params.noiseBoxHalfWidth, data=image, mask=None).data
    # if plot:
    #     imshow(noise_box)
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    logging.info("Noise statistics:\n\tMean = {:.3}\n\tMedian = {:.3}\n\tStdDev = {:.3}".format(mean, median, std))

    # Find stars
    daofind = DAOStarFinder(fwhm=params.starfinderFwhm, threshold=params.noiseThreshold*std)
    if background_subtraction:
        logging.info("Subtracting background...")
        image -= median
    logging.info("Finding sources...")
    sources = daofind(image)
    sources.sort('peak', reverse=True)
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
    logging.info("Found {} sources:\n{}".format(len(sources), sources))


    # Write results to file
    logging.info("Writing list of sources to file {}".format(params.allStarsFile))
    sources.write(params.allStarsFile, format='ascii', overwrite=True)

    # Plot results
    if plot:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.)
        norm = ImageNormalize(stretch=SqrtStretch())
        plt.imshow(image, cmap='Greys', origin='lower', norm=norm)
        plt.colorbar(pad=0.0)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.show()
        plt.close()


if __name__ == '__main__':
    main()
