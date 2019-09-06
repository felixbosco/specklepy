import argparse
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import photutils as pu


try:
    from holopy.io.paramhandler import ParamHandler
    from holopy.core.aperture import Aperture
except ModuleNotFoundError:
    # Prepare import with hardcoded path
    import warnings
    PATH = '/home/bosco/Documents/phd/sowat/pipeline/github_holopy'
    warnings.warn("Importing holopy from hardcoded path {}. Apparently holopy is not installed properly on your machine!".format(PATH), ImportWarning)
    import sys
    sys.path.insert(0, PATH)

    # Repeat import
    from holopy.io.paramhandler import ParamHandler
    from holopy.core.aperture import Aperture


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
    defaults_file = "holopy/config/holography_defaults.cfg"
    essential_attributes = ['allStarsFile', 'noiseBoxX', 'noiseBoxY', 'noiseBoxHalfWidth', 'noiseThreshold']

    if args.file is None:
        logging.error("No file or file list was provided! Use --help for instructions.")
        raise RuntimeError("No file or file list was provided! Use --help for instructions.")

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    else:
        params = ParamHandler(parameter_file=args.parameter_file, defaults_file=defaults_file, essential_attributes=essential_attributes)

    image = fits.getdata(args.file)
    noise_box = Aperture(params.noiseBoxX, params.noiseBoxY, params.noiseBoxHalfWidth, data=image)
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    print(mean, median, std)


if __name__ == '__main__':
    main()
