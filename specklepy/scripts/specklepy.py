#!/usr/bin/env python

import os

from specklepy.core.holography import holography
from specklepy.core.ssa import ssa
from specklepy.io.argparser import GeneralArgParser
from specklepy.io.parameterset import ReductionParameterSet, HolographyParameterSet
from specklepy.logging import logger
from specklepy.reduction import setup, run
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects


def main():

    # Parse args
    parser = GeneralArgParser()
    args = parser.parse_args()

    if args.debug:
        logger.setLevel('DEBUG')
        logger.debug(args)

    # Execute the script of the corresponding command
    if args.command is 'generate':

        # Read parameters from file and generate exposures
        objects = get_objects(args.parfile, debug=args.debug)
        kwargs = objects['parameters']
        generate_exposure(target=objects['target'],
                          telescope=objects['telescope'],
                          detector=objects['detector'],
                          debug=args.debug,
                          **kwargs)

    elif args.command is 'reduce':

        # In setup mode
        if args.setup:
            setup.setup(path=args.path, instrument=args.instrument, parfile=args.parfile,
                        filelist=args.filelist, sortby=args.sortby)
            return 0  # Quit program for interaction with the new parameter file

        # Else start reduction following the parameter file
        params = ReductionParameterSet(args.parfile)
        run.all(params, debug=args.debug)

    elif args.command is 'ssa':

        # Prepare path information and execute reconstruction
        if args.tmpdir is not None and not os.path.isdir(args.tmpdir):
            os.mkdir(args.tmpdir)
        ssa(args.files, mode=args.mode, tmp_dir=args.tmpdir, outfile=args.outfile, debug=args.debug)

    elif args.command is 'holography':

        # Read parameters from file and execute reconstruction
        params = HolographyParameterSet(args.parfile)
        holography(params, mode=params.options.reconstructionMode, debug=args.debug)

    elif args.command is 'aperture':
        pass
