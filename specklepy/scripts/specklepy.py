#!/usr/bin/env python

import os

from specklepy.core.holography import holography
from specklepy.core.ssa import ssa
from specklepy.io.argparser import GeneralArgParser
from specklepy.io import config
from specklepy.io.parameterset import ReductionParameterSet, HolographyParameterSet
from specklepy.logging import logger
from specklepy.reduction import setup, run
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects
from specklepy.gui.window import start


def main():

    # Parse args
    parser = GeneralArgParser()
    args = parser.parse_args()

    if args.debug:
        logger.setLevel('DEBUG')
        logger.debug(args)

    if args.gui:
        start()

    # Execute the script of the corresponding command
    if args.command is 'generate':

        # Read parameters from file and generate exposures
        target, telescope, detector, kwargs = get_objects(args.parfile, debug=args.debug)
        generate_exposure(target=target, telescope=telescope, detector=detector, debug=args.debug, **kwargs)

    elif args.command is 'reduce':

        # In setup mode
        if args.setup:
            run.setup(path=args.path, instrument=args.instrument, par_file=args.parfile, list_file=args.filelist,
                      sort_by=args.sortby)
        # Else start reduction following the parameter file
        else:
            params = config.read(args.parfile)
            run.full_reduction(params, debug=args.debug)

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
