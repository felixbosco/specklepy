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
                        '("speckle") data. The program is able to generate mock images, '
                        'reduce data, reconstruct images and to analyse the final images. Two '
                        'algorithms are implemented for image reconstruction: the SSA and '
                        'Holography algorithm.',
            epilog="Execute 'specklepy <command> -h' for further information on the commands.")
        # self.parser.add_argument('--gui', action='store_true', help='Start the graphical user interface (GUI).')
        self.parser.set_defaults(gui=False)
        self.parser.add_argument('-d', '--debug', action='store_true', help='Show debugging information.')
        subparsers = self.parser.add_subparsers(help='Available commands in Specklepy:')

        # Parser for data reduction
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
        parser_reduction.add_argument('-r', '--recursive', action='store_true',
                                      help="Search for files recursively. Note that you should provide the path in "
                                           "quotation marks and the '/**/' marker to indicate, where the recursion "
                                           "should start.")
        parser_reduction.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for SSA reconstruction
        parser_ssa = subparsers.add_parser('ssa', help='Image reconstruction with the SSA algorithm.')
        parser_ssa.set_defaults(command='ssa')
        parser_ssa.add_argument('files', nargs='+',
                                help='Generic FITS file name or list of files to consider for the SSA reconstruction.')
        parser_ssa.add_argument('-o', '--outfile', type=str, default='ssa.fits', help='Name of the output file.')
        parser_ssa.add_argument('-m', '--mode', type=str, default='same', choices=['same', 'full', 'valid'],
                                help="The mode defines the final image size. In 'same' mode, the final image will have "
                                     "the size of the first input image. In 'full' mode, every patch of the sky that "
                                     "is covered by at least one exposure will be part of the output image. In 'valid' "
                                     "mode, the output image will only cover the cross section of all exposures.")
        parser_ssa.add_argument('-r', '--reference', type=str, default=None,
                                help="Name of the reference file, serving as reference frame for frame alignment.")
        parser_ssa.add_argument('-a', '--aperture_radius', type=int, default=None,
                                help='Radius of an aperture for constraining the search of the emission peak during '
                                     'frame alignment. The center of the aperture is estimated based on a graphical '
                                     'selection of a reference star.')
        parser_ssa.add_argument('-c', '--collapse', dest='integration_method', action='store_const', const='collapse',
                                default='ssa',
                                help='Collapse the individual data cubes instead of using SSA. This option is useful '
                                     'for very faint objects.')
        parser_ssa.add_argument('-s', '--sources', dest='alignment_method', action='store_const', const='sources',
                                default='correlation',
                                help='Use sources in the field instead of a complete correlation. This is advantageous '
                                     'in the presence of hot pixels or similar features.')
        parser_ssa.add_argument('--mask', action='store_true', help='Mask hot pixels prior to aligning frames.')
        parser_ssa.add_argument('--mask_file', type=str, help='Name of a file containing a custom bad pixel mask.')
        parser_ssa.add_argument('-t', '--tmpdir', type=str, default='tmp/', help='Path for saving temporary files.')
        parser_ssa.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for holographic reconstruction
        parser_holography = subparsers.add_parser('holography',
                                                  help='Image reconstruction with the Holography algorithm.')
        parser_holography.set_defaults(command='holography')
        parser_holography.add_argument('parfile', type=str, help='Path to a parameter file.')
        parser_holography.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for aperture analysis
        parser_aperture = subparsers.add_parser('aperture', help='Aperture analysis in the image data.')
        parser_aperture.set_defaults(command='aperture')
        # parser_aperture.add_argument('mode', choices=['psf1d', 'variance', 'all'], default='all',
        #                              help='Modes for the aperture analysis')
        parser_aperture.add_argument('file', type=str, default=None, help='File name to analyse.')
        parser_aperture.add_argument('-i', '--index', nargs='+', type=int,
                                     help='Center index of the aperture to analyse. Parse this as "--index 123 456". '
                                          'Note that the order of indexes is y-x!')
        parser_aperture.add_argument('-r', '--radius', type=int, default=10,
                                     help='Radius of the aperture to analyse in pix.')
        parser_aperture.add_argument('-p', '--pixel_scale', type=float, default=1,
                                     help='Pixel scale of the data in units of arcsec.')
        parser_aperture.add_argument('-c', '--centering', action='store_true',
                                     help='Centers the aperture on the intensity peak.')
        parser_aperture.add_argument('-n', '--normalize', type=str, default=None,
                                     help='Normalize the flux values, to either "peak", "aperture" or leave as `None` '
                                          'for not normalizing.')
        parser_aperture.add_argument('-o', '--out_file', type=str, default=None, help='Saves the results to this file.')
        parser_aperture.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for plotting
        parser_plot = subparsers.add_parser('plot', help='Plot data from an input file.')
        parser_plot.set_defaults(command='plot')
        parser_plot.add_argument('file', type=str, default=None, help='Name of a file with image or table data.')
        parser_plot.add_argument('-e', '--extension', type=str, default=None, help='Name of a FITS file extension.')
        parser_plot.add_argument('-f', '--format', type=str, default=None, help='Format of the table.')
        parser_plot.add_argument('-c', '--columns', type=str, default=None, nargs=2, help='Columns of the table.')
        parser_plot.add_argument('-l', '--layout', type=str, default=None,
                                 help="Apply a default paper layout. Options are 'text' and 'column' for the requested "
                                      "image width.")
        parser_plot.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for inspecting FITS headers
        parser_inspect = subparsers.add_parser('inspect', help='Inspect FITS headers and list attributes.')
        parser_inspect.set_defaults(command='inspect')
        parser_inspect.add_argument('files', type=str, default=None, nargs='+', help='List of FITS file names.')
        parser_inspect.add_argument('-k', '--keywords', type=str, default=None, nargs='+',
                                    help='List of FITS header keywords.')
        parser_inspect.add_argument('-s', '--save', type=str, default=None,
                                    help='Name of the output file containing the table.')
        parser_inspect.add_argument('-r', '--recursive', action='store_true', help='Search for files recursively.')
        parser_inspect.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for differentiating FITS cubes
        parser_diff = subparsers.add_parser('diff', help='Differentiate cubes along time axis, for post-correlation.')
        parser_diff.set_defaults(command='diff')
        parser_diff.add_argument('files', type=str, default=None, nargs='+', help='List of FITS file names.')
        parser_diff.add_argument('-d', '--delta', type=int, default=1,
                                 help='Difference in frame number of differentiated frames. By default, every frame is '
                                      'subtracted from the subsequent one. For creating longer exposures, however, you '
                                      'can also set delta to larger values and subtract only every n-th frame. The '
                                      'resulting data cube will be correspondingly smaller.')
        parser_diff.add_argument('-k', '--keyword', type=str, default=None,
                                 help='Common part of header keywords for the individual exposure time stamp.')
        parser_diff.add_argument('-l', '--linear_regression', action='store_true',
                                 help='Apply linear regression instead of straight forward differentiation.')
        parser_diff.add_argument('-e', '--extension', type=str, default=None, help='Extension of the FITS file.')
        parser_diff.add_argument('--dtype', type=str, default=None, help='Data type to cast the data to before diff.')
        parser_diff.add_argument('--debug', action='store_true', help='show debugging information.')

        # Parser for quick apodization computations
        parser_apodization = subparsers.add_parser('apodization', help='Compute the apodization function parameters.')
        parser_apodization.set_defaults(command='apodization')
        parser_apodization.add_argument('diameter', type=float,
                                        help='diameter of the telescope primary mirror in units of meters')
        parser_apodization.add_argument('wavelength', type=float,
                                        help='observing wavelength in units of meters')
        parser_apodization.add_argument('-p', '--pixel_scale', type=float, default=None,
                                        help='pixel/ plate scale of the detector in units of mas')
        parser_apodization.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for source extraction
        parser_extraction = subparsers.add_parser('extract', help='Extract sources from an image.')
        parser_extraction.set_defaults(command='extract')
        parser_extraction.add_argument('file_name', type=str, help='Name of the image file.')
        parser_extraction.add_argument('-n', '--noise_threshold', type=float, default=5,
                                       help='Multiple of the image uncertainty.')
        parser_extraction.add_argument('-f', '--fwhm', type=float, default=5,
                                       help='Expected source FWHM in units of pixels.')
        parser_extraction.add_argument('-a', '--algorithm', default='DAO',
                                       help="Name of the StarFinder algorithm. Can be 'IRAF', 'DAO' (default), or "
                                            "'PeakFinder'.")
        parser_extraction.add_argument('-v', '--var', type=str, default=None,
                                       help='Value of the image variance or name of the FITS file extension containing '
                                            'the variance.')
        parser_extraction.add_argument('-s', '--show', action='store_true', help='Show the identified stars.')
        parser_extraction.add_argument('-o', '--out_file', type=str, default=None,
                                       help='Name of the file to store the result in.')
        parser_extraction.add_argument('--dtype', type=str, default=None,
                                       help='Cast the data type of the input image to another type.')
        parser_extraction.add_argument('-c', '--collapse', action='store_true',
                                       help='Collapse a data cube along the third axis.')
        parser_extraction.add_argument('--select', type=str, default=False,
                                       help='Provide a file name to visually select identified stars and save them to '
                                            'this file.')
        parser_extraction.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for generating mock images
        parser_generate = subparsers.add_parser('generate', help='Generate mock exposures.')
        parser_generate.set_defaults(command='generate')
        parser_generate.add_argument('parfile', type=str, help='Path to a parameter file.')
        parser_generate.add_argument('-d', '--debug', action='store_true', help='show debugging information.')

        # Parser for generating star tables
        parser_star_table = subparsers.add_parser('startable', help='Generate a star table')
        parser_star_table.set_defaults(command='startable')
        parser_star_table.add_argument('number_stars', type=int, help='Number of stars that shall be simulated.')
        parser_star_table.add_argument('out_file', type=str, help='Name of the file to store the table in.')
        parser_star_table.add_argument('-i', '--iso_files', nargs='+', type=str, help='Name(s) of the isochrone files.')
        parser_star_table.add_argument('-l', '--lf_files', nargs='+', type=str,
                                       help='Name(s) of the luminosity function files.')
        parser_star_table.add_argument('-b', '--lf_band', type=str,
                                       help='Name of the spectral band that is represented in the luminosity files.')
        parser_star_table.add_argument('-r', '--half_light_radius', type=str,
                                       help='Half-light radius of the generated cluster in units of parsecs.')
        parser_star_table.add_argument('-w', '--population_weights', nargs='+', type=float, default=None,
                                       help='Fractions of the total number of stars that shall be allocated to a given '
                                            'stellar population, represented by an isochrone and luminosity function'
                                            'file.')
        parser_star_table.add_argument('-s', '--seed', type=int, default=None,
                                       help='Seed for the random number generator.')
        parser_star_table.add_argument('-f', '--table_format', type=str, default=None,
                                       help='Format string for formatting ASCII table files.')
        parser_star_table.add_argument('-o', '--overwrite', action='store_true',
                                       help='Allow for overwriting existing files.')

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
