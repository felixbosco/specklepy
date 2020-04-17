import argparse


class GeneralArgParser(object):

    """General argument parser.

    This argument parser defines argument parsing to the general script 'specklepy'.

    Methods:
        parse_args(*args, **kwargs):
            Parse the command line input.
    """

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Specklepy is a program for the analysis of astronomical short-exposure '
                        '("speckle") data. The program is able to generate synthetic images, '
                        'reduce data, reconstruct images and to analyse the final images. Two '
                        'algorithms are implemented for image reconstruction: the SSA and '
                        'Holography algorithm.',
            epilog="Execute 'specklepy <command> -h' for further information on the commands.")
        self.parser.add_argument('-d', '--debug', action='store_true', help='show debugging information.')
        subparsers = self.parser.add_subparsers(help='Available commands in Specklepy:')

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
        parser_holography = subparsers.add_parser('holography',
                                                  help='Image reconstruction with the Holography algorithm.')
        parser_holography.set_defaults(command='holography')
        parser_holography.add_argument('parfile', type=str, help='Path to a parameter file')

        # Parser for aperture analysis
        parser_aperture = subparsers.add_parser('aperture', help='Aperture analysis in the image data.')
        parser_aperture.set_defaults(command='aperture')
        parser_aperture.add_argument('mode', choices=[], help='Modes for the aperture analysis')


    def parse_args(self, *args, **kwargs):
        """Parse command line arguments.

        Args:
            *args:
                Parsed to self.parser.parse_args().
            **kwargs:
                Parsed to self.parser.parse_args().

        Returns:
            args (argparse.Namespace):
                argparse.Namespace object that contains parsed arguments from self.parser.parse_args().
        """
        return self.parser.parse_args(*args, **kwargs)
