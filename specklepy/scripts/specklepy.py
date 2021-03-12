#!/usr/bin/env python3

import os
import sys

from specklepy.core import analysis
from specklepy.core.holography import holography
from specklepy.core.sourceextraction import extract_sources
from specklepy.core.ssa import ssa
from specklepy.io.argparser import GeneralArgParser
from specklepy.io import config
from specklepy.logging import logger
from specklepy.plotting.plot import Plot
from specklepy.reduction import diff
from specklepy.reduction import run
from specklepy.reduction.datareduction import DataReduction
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects
from specklepy.utils.resolution import get_resolution_parameters
from specklepy.gui.window import start


def main():

    # Initial notification
    logger.debug("Starting Specklepy...")

    # Parse args
    logger.debug("Parsing terminal input...")
    parser = GeneralArgParser()
    args = parser.parse_args()
    try:
        logger.debug(f"Specklepy command: {args.command}")
    except AttributeError:
        sys.tracebacklimit = 0
        raise RuntimeError("No command was provided. Execute 'specklepy -h' for help!")
    logger.debug(f"Parsed terminal input:\n{args}")

    # Set logging level
    if args.debug:
        logger.setLevel('DEBUG')
        logger.debug(args)

    # Start the GUI
    if args.gui:
        start()

    # Execute the script of the corresponding command
    if args.command == 'aperture':
        if args.mode == 'all':
            logger.info("Extract aperture statistics")
            analysis.aperture_analysis(file=args.file, index=args.index, radius=args.radius, out_file=args.out_file,
                                       pixel_scale=args.pixel_scale, recenter=args.centering, debug=args.debug)
        elif args.mode == 'psf1d':
            logger.info("Extract 1D PSF profile")
            analysis.get_psf_1d(args.file, args.index, args.radius, args.out_file, args.normalize, debug=args.debug)
        elif args.mode == 'variance':
            logger.info("Extract 1D PSF variation")
            analysis.get_psf_variation(args.file, args.index, args.radius, args.out_file, args.normalize, args.debug)
        else:
            logger.warning(f"Aperture mode {args.mode} not recognized!")

    elif args.command == 'apodization':
        get_resolution_parameters(wavelength=args.wavelength, diameter=args.diameter, pixel_scale=args.pixel_scale)

    elif args.command == 'diff':
        # Set differentiation method
        if args.linear_regression:
            method = 'linreg'
        else:
            method = 'direct'
        diff.differentiate_cube(files=args.files, delta=args.delta, method=method, exposure_time_prefix=args.keyword,
                                extension=args.extension, dtype=args.dtype)

    elif args.command == 'extract':
        if args.out_file is None:
            args.out_file = 'default'
        extract_sources(image=args.file_name, noise_threshold=args.noise_threshold, fwhm=args.fwhm, image_var=args.var,
                        show=args.show, collapse=args.collapse, write_to=args.out_file, cast_dtype=args.dtype,
                        select=args.select, debug=args.debug)

    elif args.command == 'generate':

        # Read parameters from file and generate exposures
        target, telescope, detector, parameters = get_objects(args.parfile, debug=args.debug)
        generate_exposure(target=target, telescope=telescope, detector=detector, debug=args.debug, **parameters)

    elif args.command == 'holography':

        # Read parameters from file and execute reconstruction
        defaults_file = os.path.join(os.path.dirname(__file__), '../config/holography.cfg')
        defaults_file = os.path.abspath(defaults_file)
        params = config.read(defaults_file)
        params = config.update_from_file(params, args.parfile)
        holography(params, mode=params['OPTIONS']['reconstructionMode'], debug=args.debug)

    elif args.command == 'inspect':
        run.inspect(files=args.files, keywords=args.keywords, save=args.save, recursive=args.recursive,
                    debug=args.debug)

    elif args.command == 'plot':
        plot = Plot.from_file(file_name=args.file, extension=args.extension, columns=args.columns, format=args.format,
                              layout=args.layout, debug=args.debug)
        plot.apply_layout(layout=args.layout)
        plot.save()
        plot.show()

    elif args.command == 'reduce':

        # In setup mode
        if args.setup:
            DataReduction.set_up(path=args.path, instrument=args.instrument, par_file=args.parfile,
                                 list_file=args.filelist, sort_by=args.sortby, recursive=args.recursive)
        # Else start reduction following the parameter file
        else:
            # params = config.read(args.parfile)
            # run.full_reduction(params, debug=args.debug)
            data_reduction = DataReduction.from_file(args.parfile)
            data_reduction.run()

    elif args.command == 'ssa':

        # Prepare path information and execute reconstruction
        if args.tmpdir is not None and not os.path.isdir(args.tmpdir):
            os.mkdir(args.tmpdir)
        ssa(files=args.files, mode=args.mode, reference_file=args.reference, tmp_dir=args.tmpdir, outfile=args.outfile,
            box_indexes=args.box_indexes, integration_method=args.integration_method, debug=args.debug)
