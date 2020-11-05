#!/usr/bin/env python

import os

from specklepy.core import analysis
from specklepy.core.holography import holography
from specklepy.core.sourceextraction import extract_sources
from specklepy.core.ssa import ssa
from specklepy.io.argparser import GeneralArgParser
from specklepy.io import config
from specklepy.logging import logger
from specklepy.plotting.plot import Plot
from specklepy.reduction import run
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects
from specklepy.utils.resolution import get_resolution_parameters
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
        target, telescope, detector, parameters = get_objects(args.parfile, debug=args.debug)
        generate_exposure(target=target, telescope=telescope, detector=detector, debug=args.debug, **parameters)

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
        ssa(args.files, mode=args.mode, tmp_dir=args.tmpdir, outfile=args.outfile, box_indexes=args.box_indexes,
            debug=args.debug)

    elif args.command is 'holography':

        # Read parameters from file and execute reconstruction
        defaults_file = os.path.join(os.path.dirname(__file__), '../config/holography.cfg')
        defaults_file = os.path.abspath(defaults_file)
        params = config.read(defaults_file)
        params = config.update_from_file(params, args.parfile)
        holography(params, mode=params['OPTIONS']['reconstructionMode'], debug=args.debug)

    elif args.command is 'aperture':
        if args.mode == 'psf1d':
            logger.info("Extract 1D PSF profile")
            analysis.get_psf_1d(args.file, args.index, args.radius, args.out_file, args.normalize, debug=args.debug)
        elif args.mode == 'variance':
            logger.info("Extract 1D PSF variation")
            analysis.get_psf_variation(args.file, args.index, args.radius, args.out_file, args.normalize, args.debug)
        else:
            logger.warning(f"Aperture mode {args.mode} not recognized!")

    elif args.command is 'extract':
        if args.out_file is None:
            args.out_file = 'sources_' + args.file_name.replace('.fits', '.dat')
        extract_sources(image=args.file_name, noise_threshold=args.noise_threshold, fwhm=args.fwhm, image_var=args.var,
                        write_to=args.out_file)

    elif args.command == 'plot':
        plot = Plot.from_file(file_name=args.file, extension=args.extension, columns=args.columns, debug=args.debug)
        plot.save()
        plot.show()

    elif args.command is 'apodization':
        get_resolution_parameters(wavelength=args.wavelength, diameter=args.diameter, pixel_scale=args.pixel_scale)
