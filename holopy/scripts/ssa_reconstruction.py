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
    parser.add_argument('-t', '--tmpdir', type=str, default=None, help='Path to save temporary files to.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the outfile file.')

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

    if args.tmpdir is not None and not os.path.isdir(args.tmpdir):
        os.system('mkdir {}'.format(args.tmpdir))

    # Execute reconstruction
    outfile = Outfile(file_list=files, filename=args.outfile, cards={"RECONSTRUCTION": "SSA"})
    # algorithm = SSAReconstruction()
    # algorithm.execute(files, outfile=outfile)
    ssa(files, tmp_dir=args.tmpdir, outfile=outfile)


if __name__ == '__main__':
    main()
