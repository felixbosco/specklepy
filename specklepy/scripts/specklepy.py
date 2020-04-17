#!/usr/bin/env python

from specklepy.io.argparser import GeneralArgParser
from specklepy.reduction import steps
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects


def main():

    # Parse args
    parser = GeneralArgParser()
    args = parser.parse_args()

    # Execute the script of the corresponding command
    if args.command is 'generate':

        # Read parameters from file
        objects = get_objects(args.parfile, debug=args.debug)
        kwargs = objects['kwargs']

        # Generate exposures
        generate_exposure(target=objects['target'],
                          telescope=objects['telescope'],
                          detector=objects['detector'],
                          debug=args.debug,
                          **kwargs)

    elif args.command is 'reduce':

        # In setup mode
        if args.setup:
            steps.setup(args)

            # Quit program for interactiong with the new parameter file
            return 0

        # Else start reduction following the parameter file
        pass


    elif args.command is 'ssa':
        pass

    elif args.command is 'holography':
        pass

    elif args.command is 'aperture':
        pass
