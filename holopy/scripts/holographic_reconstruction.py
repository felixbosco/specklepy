import argparse

try:
    from holopy.io.filehandler import FileHandler
    from holopy.io.outfile import Outfile
    from holopy.logging import logging
    from holopy.core.ssa import SSAReconstructor
except ModuleNotFoundError:
    # Prepare import with hardcoded path
    import warnings
    PATH = '/home/bosco/Documents/phd/sowat/pipeline/github_holopy'
    warnings.warn("Importing spampy from hardcoded path {}. Apparently holopy is not installed properly on your machine!".format(PATH), ImportWarning)
    import sys
    sys.path.insert(0, PATH)

    # Repeat import
    from holopy.io.filehandler import FileHandler
    from holopy.io.outfile import Outfile
    from holopy.logging import logging
    from holopy.core.ssa import SSAReconstructor



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

    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")

    # Execute reconstruction
    pass


if __name__ == '__main__':
    main()
