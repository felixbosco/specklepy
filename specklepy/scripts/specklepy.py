#!/usr/bin/env python

from specklepy.io.argparser import GeneralArgParser
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects


def main():

    # Parse args
    parser = GeneralArgParser()
    args = parser.parse_args()

    # Execute the script of the corresponding command
    if args.command is 'generate':

        # Read parameters from file
        if args.parfile is None:
            raise RuntimeError("No parameter file was provided! Use --help for instructions.")
        objects = get_objects(args.parfile, debug=args.debug)

        # Generate exposures
        generate_exposure(target=objects['target'],
                          telescope=objects['telescope'],
                          detector=objects['detector'],
                          **objects['kwargs'],
                          debug=args.debug)

    elif args.command is 'reduce':
        pass

    elif args.command is 'ssa':
        pass

    elif args.command is 'holography':
        pass

    elif args.command is 'aperture':
        pass
    