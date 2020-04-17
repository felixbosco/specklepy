#!/usr/bin/env python

import argparse


# Define argument parser
parser = argparse.ArgumentParser(description='Specklepy is a program for the analysis of astronomical short-exposure '
                                             '("speckle") data. The program is able to generate synthetic images, '
                                             'reduce data, reconstruct images and to analyse the final images. Two '
                                             'algorithms are implemented for image reconstruction: the SSA and '
                                             'Holography algorithm.',
                                 epilog="Execute 'specklepy <command> -h' for further information on the commands.")
parser.add_argument('-d', '--debug', action='store_true', help='show debugging information.')
subparsers = parser.add_subparsers(help='Available commands in Specklepy:')

# Parser for generating synthetic images
parser_generate = subparsers.add_parser('generate', help='Generate synthetic exposures.')
parser_generate.set_defaults(command='generate')
parser_generate.add_argument('parfile', type=str, help='Path to a parameter file')

# Parser for reduction
parser_reduction = subparsers.add_parser('reduce', help='Data reduction.')
parser_reduction.set_defaults(command='reduce')
parser_reduction.add_argument('parfile', type=str, help='Path to a parameter file')

# Parser for SSA reconstruction
parser_ssa = subparsers.add_parser('ssa', help='Image reconstruction with the SSA algorithm.')
parser_ssa.set_defaults(command='ssa')
parser_ssa.add_argument('parfile', type=str, help='Path to a parameter file')

# Parser for SSA reconstruction
parser_holography = subparsers.add_parser('holography', help='Image reconstruction with the Holography algorithm.')
parser_holography.set_defaults(command='holography')
parser_holography.add_argument('parfile', type=str, help='Path to a parameter file')

# Parser for aperture analysis
parser_aperture = subparsers.add_parser('aperture', help='Aperture analysis in the image data.')
parser_aperture.set_defaults(command='aperture')
parser_aperture.add_argument('mode', choices=[], help='Modes for the aperture analysis')


def main():

    # Parse args
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