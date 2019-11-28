#!/usr/bin/env python

import argparse
import os
import sys

try:
    from holopy.io.filehandler import FileHandler
    from holopy.io.outfile import Outfile
    from holopy.logging import logging
    # from holopy.algorithms.ssa import SSAReconstruction
except ModuleNotFoundError:
    # Prepare import with hardcoded path
    import warnings
    PATH = os.getcwd()
    warnings.warn("Importing spampy from hardcoded path {}. Apparently holopy is not installed properly on your machine!".format(PATH), ImportWarning)
    sys.path.insert(0, PATH)

    # Repeat import
    from holopy.io.filehandler import FileHandler
    from holopy.io.outfile import Outfile
    from holopy.logging import logging
    # from holopy.algorithms.ssa import SSAReconstruction
    from holopy.core.ssa import ssa



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script creates a simple shift-and-add (SSA) reconstruction of the input files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='Fits file or generic file name to consider for the SSA reconstruction.')
    parser.add_argument('-p', '--parameter_file', type=str, help='Path to the parameter file.')
    parser.add_argument('-o', '--output', type=str, default=None, help='Name of the output file.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(options=None):

    args = parser(options=options)

    if args.file is not None:
        files = FileHandler(args.file)()
    else:
        logging.error("No file or file list was provided! Use --help for instructions.")
        raise RuntimeError("No file or file list was provided! Use --help for instructions.")

    if args.parameter_file is None:
        logging.warning("No parameter file was provided! Reconstruction is executed with default values.")

    # Execute reconstruction
    outfile = Outfile(file_list=files, filename=args.output, cards={"RECONSTRUCTION": "SSA"})
    # algorithm = SSAReconstruction()
    # algorithm.execute(files, outfile=outfile)
    ssa(files, outfile=outfile)


if __name__ == '__main__':
    main()
