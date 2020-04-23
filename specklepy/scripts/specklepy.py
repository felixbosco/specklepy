#!/usr/bin/env python

# TODO: Export all the steps into external functions to make this script short and concise

import os
from astropy.utils import iers  # Suppress downloading from IERS
iers.conf.auto_download = False  # Suppress downloading from IERS

from specklepy.core.holography import holography
from specklepy.core.ssa import ssa
from specklepy.io.argparser import GeneralArgParser
from specklepy.io.filemanager import FileManager
from specklepy.io.parameterset import ParameterSet
from specklepy.logging import logger
from specklepy.reduction import setup, run
from specklepy.synthetic.generate_exposure import generate_exposure, get_objects


def main():

    # Parse args
    parser = GeneralArgParser()
    args = parser.parse_args()

    if args.debug:
        logger.debug(args)

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

        # TODO: Export reduction steps into functions to shorten this section here.

        # In setup mode
        if args.setup:
            setup.setup(path=args.path, instrument=args.instrument, parfile=args.parfile,
                        filelist=args.filelist, sortby=args.sortby)

            # Quit program for interaction with the new parameter file
            return 0

        # Else start reduction following the parameter file
        # Default values
        defaults_file = os.path.join(os.path.dirname(__file__), '../config/reduction.cfg')
        essential_attributes = {'paths': ['filePath', 'fileList', 'outDir', 'tmpDir', 'filePrefix'],
                                'setup': ['setupKeywords'],
                                'flat': ['skip', 'masterFlatFile'],
                                'sky': ['skip', 'source', 'method', 'sigmaClip']}

        # Read parameters for parameter file
        params = ParameterSet(parameter_file=args.parfile,
                              defaults_file=defaults_file,
                              essential_attributes=essential_attributes,
                              make_dirs=['tmpDir'])

        run.all(params)


    elif args.command is 'ssa':

        # Prepare path information
        files = FileManager(args.files).files
        if args.tmpdir is not None and not os.path.isdir(args.tmpdir):
            os.mkdir(args.tmpdir)

        # Execute reconstruction
        ssa(files, mode=args.mode, tmp_dir=args.tmpdir, outfile=args.outfile, debug=args.debug)

    elif args.command is 'holography':

        # Default values
        defaults_file = os.path.join(os.path.dirname(__file__), '../config/holography.cfg')
        essential_attributes = {'paths': ['inDir', 'tmpDir', 'outFile', 'alignmentReferenceFile', 'refSourceFile'],
                                'starfinder': ['starfinderFwhm', 'noiseThreshold'],
                                'psfextraction': ['mode', 'psfRadius', 'noiseThreshold', 'noiseReferenceMargin'],
                                'apodization': ['apodizationWidth', 'apodizationType'],
                                'options': ['reconstructionMode', 'varianceExtensionName']}
        make_dirs = ['tmpDir']

        # Read parameters from file
        params = ParameterSet(parameter_file=args.parfile, defaults_file=defaults_file,
                              essential_attributes=essential_attributes, make_dirs=make_dirs)

        # Execute reconstruction
        holography(params, mode=params.options.reconstructionMode, debug=args.debug)

    elif args.command is 'aperture':
        pass
