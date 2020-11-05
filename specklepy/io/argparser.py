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
        # self.parser.add_argument('--gui', action='store_true', help='Start the graphical user interface (GUI).')
        self.parser.set_defaults(gui=False)
        self.parser.add_argument('-d', '--debug', action='store_true', help='Show debugging information.')
        subparsers = self.parser.add_subparsers(help='Available commands in Specklepy:')

        # Parser for generating synthetic images
        parser_generate = subparsers.add_parser('generate', help='Generate synthetic exposures.')
        parser_generate.set_defaults(command='generate')
        parser_generate.add_argument('parfile', type=str, help='Path to a parameter file.')
        parser_generate.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for reduction
        parser_reduction = subparsers.add_parser('reduce', help='Data reduction.')
        parser_reduction.set_defaults(command='reduce')
        parser_reduction.add_argument('parfile', type=str,
                                      help='Path to a parameter file. This will be created if in setup mode')
        parser_reduction.add_argument('--setup', action='store_true', help='Switch to setup mode.')
        parser_reduction.add_argument('-p', '--path', type=str, default=None,
                                      help='Path to the files that will be listed in the outfile (in setup mode).')
        parser_reduction.add_argument('-i', '--instrument', type=str,
                                      help='Name of the instrument (for setup mode only).')
        parser_reduction.add_argument('-f', '--filelist', type=str, default='specklepy_reduction_files.tab',
                                      help="Name of the file containing the file names (for setup mode only).")
        parser_reduction.add_argument('-s', '--sortby', type=str, default=None,
                                      help="Header card to sort the output table by (for setup mode only).")
        parser_reduction.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for SSA reconstruction
        parser_ssa = subparsers.add_parser('ssa', help='Image reconstruction with the SSA algorithm.')
        parser_ssa.set_defaults(command='ssa')
        parser_ssa.add_argument('files', nargs='+',
                                help='Generic FITS file name or list of files to consider for the SSA reconstruction.')
        parser_ssa.add_argument('-m', '--mode', type=str, default='same', choices=['same', 'full', 'valid'],
                                help="The mode defines the final image size. In 'same' mode, the final image will have "
                                     "the size of the first input image. In 'full' mode, every patch of the sky that "
                                     "is covered by at least one exposure will be part of the output image. In 'valid' "
                                     "mode, the output image will only cover the cross section of all exposures.")
        parser_ssa.add_argument('-b', '--box_indexes', type=int, nargs=4, default=None,
                                help='Coordinates of a box constraining the search of the emission peak for frame '
                                     'alignment. Provide as a list [x_min, x_max, y_min, y_max].')
        parser_ssa.add_argument('-o', '--outfile', type=str, default='ssa.fits', help='Name of the output file.')
        parser_ssa.add_argument('-t', '--tmpdir', type=str, default='tmp/', help='Path for saving temporary files.')
        parser_ssa.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for holographic reconstruction
        parser_holography = subparsers.add_parser('holography',
                                                  help='Image reconstruction with the Holography algorithm.')
        parser_holography.set_defaults(command='holography')
        parser_holography.add_argument('parfile', type=str, help='Path to a parameter file.')
        parser_holography.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for extracting 1D PSF profiles
        # parser_psf_1d = subparsers.add_parser('psf1d', help='Extract 1D PSF profiles.')
        # parser_psf_1d.set_defaults(command='psf1d')
        # parser_psf_1d.add_argument('file', type=str, default=None, help='File name to analyse.')
        # parser_psf_1d.add_argument('-i', '--index', nargs='+', type=int,
        #                            help='Center index of the aperture to analyse. Provide this as "--index 123 456".')
        # parser_psf_1d.add_argument('-r', '--radius', type=int, default=10,
        #                            help='Radius of the aperture to analyse in pix.')
        # parser_psf_1d.add_argument('-n', '--normalize', type=str, default=None,
        #                            help='Normalize the flux values, to either "peak", "aperture" or leave as `None` for '
        #                                 'not normalizing.')
        # parser_psf_1d.add_argument('-o', '--out_file', type=str, default=None, help='Saves the results to this file.')
        # parser_psf_1d.add_argument('-d', '--debug', action='store_true', help='Set to inspect intermediate results.')

        # Parser for aperture analysis
        parser_aperture = subparsers.add_parser('aperture', help='Aperture analysis in the image data.')
        parser_aperture.set_defaults(command='aperture')
        parser_aperture.add_argument('mode', choices=['psf1d', 'variance'], help='Modes for the aperture analysis')
        parser_aperture.add_argument('file', type=str, default=None, help='File name to analyse.')
        parser_aperture.add_argument('-i', '--index', nargs='+', type=int,
                                   help='Center index of the aperture to analyse. Provide this as "--index 123 456".')
        parser_aperture.add_argument('-r', '--radius', type=int, default=10,
                                   help='Radius of the aperture to analyse in pix.')
        parser_aperture.add_argument('-n', '--normalize', type=str, default=None,
                                   help='Normalize the flux values, to either "peak", "aperture" or leave as `None` for '
                                        'not normalizing.')
        parser_aperture.add_argument('-o', '--out_file', type=str, default=None, help='Saves the results to this file.')
        parser_aperture.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        parser_plot = subparsers.add_parser('plot', help='Plot data from an input file.')
        parser_plot.set_defaults(command='plot')
        parser_plot.add_argument('file', type=str, default=None, help='Name of a file with image or table data.')
        parser_plot.add_argument('-e', '--extension', type=str, default=None, help='Name of a FITS file extension.')
        parser_plot.add_argument('-f', '--format', type=str, default=None, help='Format of the table.')
        parser_plot.add_argument('-c', '--columns', type=str, default=None, nargs=2, help='Columns of the table.')
        parser_plot.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for
        parser_apodization = subparsers.add_parser('apodization', help='Compute the apodization function parameters.')
        parser_apodization.set_defaults(command='apodization')
        parser_apodization.add_argument('diameter', type=float,
                                        help='diameter of the telescope primary mirror in units of meters')
        parser_apodization.add_argument('wavelength', type=float,
                                        help='observing wavelength in units of meters')
        parser_apodization.add_argument('-p', '--pixel_scale', type=float, default=None,
                                        help='pixel/ plate scale of the detector in units of mas')
        parser_apodization.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        parser_extraction = subparsers.add_parser('extract', help='Extract sources from an image.')
        parser_extraction.set_defaults(command='extract')
        parser_extraction.add_argument('file_name', type=str, help='Name of the image file.')
        parser_extraction.add_argument('-n', '--noise_threshold', type=float, help='Multiple of the image uncertainty.')
        parser_extraction.add_argument('-f', '--fwhm', type=float, help='Expected source FWHM in units of pixels.')
        parser_extraction.add_argument('-v', '--var', type=str, default=None,
                                       help='Value of the image variance or name of the FITS file extension containing '
                                            'the variance.')
        parser_extraction.add_argument('-o', '--out_file', type=str, default=None,
                                       help='Name of the file to store the result in.')

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
