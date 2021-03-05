#!/usr/bin/env python

import argparse
import os
import sys

from specklepy.logging import logger
from specklepy.io.parameterset import ParameterSet
from specklepy.io.filearchive import FileArchive
from specklepy.reduction.flat import MasterFlat
from specklepy.reduction import sky



def parser(options=None):

    parser = argparse.ArgumentParser(description='This script reduces the data, following the parameters specified in the paramater files.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('parameter_file', type=str, help='Path to the parameter file.')
    parser.add_argument('-d', '--debug', action='store_true', help='Set to inspect intermediate results.')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args



def main(options=None):

    args = parser(options=options)

    # Default values
    defaults_file = os.path.join(os.path.dirname(__file__), '../config/reduction.cfg')
    essential_attributes = {'paths': ['filePath', 'fileList', 'outDir', 'tmpDir', 'filePrefix'],
                            'setup': ['setupKeywords'],
                            'flat': ['skipFlat', 'flatCorrectionPrefix'],
                            'sky': ['skipSky', 'ignoreTimeStamps', 'skySubtractionPrefix', 'fromImage', 'backgroundMethod', 'backgroundSigmaClip']}
    make_dirs = ['tmpDir']

    # Read parameters from file
    if args.parameter_file is None:
        raise RuntimeError("No parameter file was provided! Use --help for instructions.")
    else:
        params = ParameterSet(parameter_file=args.parameter_file,
                        defaults_file=defaults_file,
                        essential_attributes=essential_attributes,
                        make_dirs=make_dirs)

    # Check parameter values
    if params.sky.backgroundMethod not in ['scalar']:
        raise RuntimeError(f"Background subtraction method {params.sky.backgroundMethod!r} is not implemented.")

    # Execute data reduction
    # (0) Read file list table
    logger.info("Reading file list ...")
    inFiles = FileArchive(params.paths.fileList)
    # inFiles.table['PATH'] = [params.paths.filePath] * len(inFiles.table['FILE'])
    # inFiles.table['CORRECTIONS'] = [''] * len(inFiles.table['FILE'])
    print(inFiles.table)
    # last_correction = None

    # (0) Make reduction files



    # (1) Flat fielding
    if not params.flat.skipFlat:
        flat_files = inFiles.filter({'OBSTYPE': 'FLAT'})
        master_flat = MasterFlat(flat_files, file_name=params.flat.masterFlatFile, file_path=params.paths.filePath)
        master_flat.combine()
        to_be_flatfield_corrected = inFiles.filter({'OBSTYPE': ['SCIENCE', 'SKY']})
        flatfield_corrected = master_flat.run_correction(to_be_flatfield_corrected, prefix=params.flat.flatCorrectionPrefix)#, savedir=params.paths.tmpDir)
        inFiles.update_filenames(flatfield_corrected)


    # (...) Linearisation

    # (...) Identify setups
    inFiles.identify_setups(params.setup.setupKeywords)

    # (...) Sky subtraction
    if not params.sky.skipSky:
        if hasattr(params.sky, 'fromImage') and params.sky.fromImage:
            logger.info("Estimating background flux from the images...")
            science_files = inFiles.filter({'OBSTYPE': 'Science'})
            if params.sky.backgroundMethod is 'scalar':
                sky.subtract_scalar_background(science_files, params=params, prefix=params.sky.skySubtractionPrefix, debug=args.debug)
            elif params.sky.backgroundMethod is '2D':
                # Todo implement subtraction of 2D background
                pass
        else:
            logger.info("Estimating background flux from sky frames...")
            sequences = sky.identify_sequences(inFiles.table, file_path=params.paths.filePath, ignore_time_stamps=params.sky.ignoreTimeStamps)
            for sequence in sequences:
                sequence.subtract_master_sky(save_to=params.paths.tmpDir, filename_prefix=params.sky.skySubtractionPrefix)



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logger.info('Interrupted by user...')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
