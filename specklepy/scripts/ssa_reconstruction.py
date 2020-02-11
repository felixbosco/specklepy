#!/usr/bin/env python

import argparse
import os
import sys

from specklepy.logging import logging
from specklepy.io.filemanager import FileManager
from specklepy.io.reconstructionfile import ReconstructionFile
from specklepy.core.ssa import ssa



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script creates a simple shift-and-add (SSA) reconstruction of the input files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-f', '--file', type=str, default=None, help='Fits file or generic file name to consider for the SSA reconstruction.')
    parser.add_argument('-t', '--tmpdir', type=str, default=None, help='Path to save temporary files to.')
    parser.add_argument('-o', '--outfile', type=str, default=None, help='Name of the output file.')
    parser.add_argument('-d', '--debug', type=bool, default=False, help='Set to True to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)

    if args.file is not None:
        files = FileManager(args.file)()
    else:
        logging.error("No file or file list was provided! Use --help for instructions.")
        raise RuntimeError("No file or file list was provided! Use --help for instructions.")

    if args.tmpdir is not None and not os.path.isdir(args.tmpdir):
        os.system('mkdir {}'.format(args.tmpdir))

    # Execute reconstruction
    ssa(files, tmp_dir=args.tmpdir, outfile=args.outfile, debug=args.debug)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logging.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
