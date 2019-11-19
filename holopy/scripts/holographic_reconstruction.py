import argparse
import os
import sys

try:
    from holopy.io.filehandler import FileHandler
    from holopy.io.outfile import Outfile
    from holopy.io.paramhandler import ParamHandler
    from holopy.logging import logging
    from holopy.algorithms.holography import HolographicReconstruction
except ModuleNotFoundError:
    # Prepare import with hardcoded path
    import warnings
    PATH = os.getcwd()
    warnings.warn("Importing holopy from hardcoded path {}. Apparently holopy is not installed properly on your machine!".format(PATH), ImportWarning)
    sys.path.insert(0, PATH)

    # Repeat import
    from holopy.io.filehandler import FileHandler
    from holopy.io.outfile import Outfile
    from holopy.io.paramhandler import ParamHandler
    from holopy.logging import logging
    from holopy.algorithms.holography import HolographicReconstruction



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script creates a simple shift-and-add (SSA) reconstruction of the input files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--parameter_file', type=str, help='Path to the parameter file.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(options=None):

    args = parser(options=options)

    # Default values
    defaults_file = "holopy/config/holography_defaults.cfg"
    essential_attributes = ['inDir', 'tmpDir', 'outFile', 'refSourceFile', 'psfRadius', 'maskRadius', 'noiseThreshold', 'apodizationWidth', 'apodizationType']
    make_dirs = ['inDir', 'tmpDir']

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    else:
        params = ParamHandler(parameter_file=args.parameter_file,
                        defaults_file=defaults_file,
                        essential_attributes=essential_attributes,
                        make_dirs=make_dirs)

    # Instantiate handler classes
    params.outFile = Outfile(file_list=params.inFiles, filename=params.outFile, cards={"RECONSTRUCTION": "Holography"})

    # Execute reconstruction
    algorithm = HolographicReconstruction(params)
    algorithm.execute()


if __name__ == '__main__':
    main()
