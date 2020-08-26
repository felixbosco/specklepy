from datetime import datetime
import glob
import os

from astropy.io import fits
from astropy.table import Table

from specklepy.io import config
from specklepy.io.filearchive import FileArchive
from specklepy.logging import logger
from specklepy.reduction import flat, sky


def setup(path, instrument, par_file=None, list_file=None, sort_by=None):
    """Sets up the data reduction parameter file and file list.

    Args:
        path (str):
            Path to the files.
        instrument (str):
            Name of the instrument that took the data. This must be covered by config/instruments.cfg.
        par_file (str, optional):
            Name of the output default parameter file for the reduction.
        list_file (str):
            Name of the output file that contains all the file names and header information.
        sort_by (str, optional):
            Header card that is used for the sorting of files.
    """

    # Defaults
    header_cards = ['OBSTYPE', 'OBJECT', 'FILTER', 'EXPTIME', 'nFRAMES', 'DATE']
    dtypes = [str, str, str, float, int, str]
    instrument_config_file = os.path.join(os.path.dirname(__file__), '../config/instruments.cfg')

    # Read config
    configs = config.read(instrument_config_file)
    instrument = configs['INSTRUMENTS'][instrument]
    instrument_header_cards = configs[instrument]

    # Double check whether all aliases are defined
    cards = []
    for card in header_cards:
        try:
            cards.append(instrument_header_cards[card])
        except KeyError:
            logger.info(
                f"Dropping header card {card} from setup identification, as there is no description in the config file."
                f"\nCheck out {instrument_config_file} for details.")
            header_cards.remove(card)

    # Apply fall back values
    if path is None:
        path = '.'
    if list_file is None:
        list_file = 'files.tab'
    if par_file is None:
        par_file = 'reduction.yaml'

    # Find files
    if '*' in path:
        files = glob.glob(path)
    else:
        files = glob.glob(os.path.join(path, '*fits'))
    if len(files):
        logger.info(f"Found {len(files)} file(s)")
        files.sort()
    else:
        logger.error(f"Found no files in {path}!")
        raise RuntimeError(f"Found no files in {path}!")

    # Initialize a file archive
    raw_files = FileArchive(files, cards=cards, dtypes=dtypes, names=header_cards)
    raw_files.identify_setups(['FILTER', 'EXPTIME'])
    # table = raw_files.table

    # # Initialize output file information table
    # table = Table(names=['FILE'] + header_cards, dtype=[str] + dtypes)
    #
    # # Read data from files
    # for file in files:
    #     logger.info(f"Retrieving header information from file {file}")
    #     hdr = fits.getheader(file)
    #     new_row = [os.path.basename(file)]
    #     for card in header_cards:
    #         try:
    #             new_row.append(hdr[instrument_header_cards[card]])
    #         except KeyError:
    #             logger.info(f"Skipping file {os.path.basename(file)} due to at least one missing header card "
    #                         f"({instrument_header_cards[card]}).")
    #             break
    #     if len(new_row) == len(table.columns):
    #         table.add_row(new_row)
    #
    # # Sort table entries by default properties and user request
    # table.sort('FILE')
    # table.sort('OBSTYPE')
    # if sort_by:
    #     table.sort(sort_by)
    #
    # # Identify instrument setups
    # setups = identify_instrument_setups(table)
    # table.add_column(setups)
    #
    # # Save table
    # logger.info(f"Writing header information to {list_file}")
    # table.write(list_file, format='ascii.fixed_width', overwrite=True)
    raw_files.write_table(file_name=list_file)

    # Write dummy parameter file for the reduction
    _, ext = os.path.splitext(par_file)
    if 'yaml' in ext:
        logger.info(f"Creating default reduction YAML parameter file {par_file}")
        par_file_content = f"PATHS:\n  filePath: {path}\n  fileList: {list_file}\n  outDir: Science/\n  tmpDir: tmp/" \
                           f"\n\nFLAT:\n  masterFlatFile: MasterFlat.fits" \
                           f"\n\nSKY:\n  method: scalar"
    else:
        logger.info(f"Creating default reduction INI parameter file {par_file}")
        par_file_content = f"[PATHS]\nfilePath = {path}\nfileList = {list_file}\noutDir = Science/\ntmpDir = tmp/" \
                           f"\n\n[FLAT]\nmasterFlatFile = MasterFlat.fits" \
                           f"\n\n[SKY]\nmethod = scalar"
    with open(par_file, 'w+') as par_file:
        par_file.write(par_file_content)


def full_reduction(params, debug=False):
    """Execute a full reduction following the parameters in the `params` dictionary.

    TODO: Split this function into the parts and sort into the other modules

    Args:
        params (dict):
            Dictionary with all the settings for reduction
        debug (bool, optional):
            Show debugging information
    """

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')

    # (0) Read file list table
    logger.info("Reading file list ...")
    in_files = FileArchive(params['PATHS']['fileList'],
                           in_dir=params['PATHS']['filePath'],
                           out_dir=params['PATHS']['outDir'])
    logger.info('\n' + str(in_files.table))

    # (1) Initialize reduction files
    if not os.path.isdir(params['PATHS']['outDir']):
        os.makedirs(params['PATHS']['outDir'])
    product_files = in_files.initialize_product_files()

    # (2) Flat fielding
    if 'skip' in params['FLAT'] and params['FLAT']['skip']:
        logger.info('Skipping flat fielding as requested from parameter file...')
    else:
        flat_files = in_files.get_flats()
        if len(flat_files) == 0:
            logger.warning("Did not find any flat field observations. No flat field correction will be applied!")
        else:
            logger.info("Starting flat field correction...")
            master_flat = flat.MasterFlat(flat_files, file_name=params['FLAT']['masterFlatFile'],
                                          file_path=params['PATHS']['filePath'])
            master_flat.combine()
            master_flat.run_correction(file_list=product_files, file_path=None)

    # (3) Linearization
    # TODO: Implement linearization

    # (4) Sky subtraction
    if 'skip' in params['SKY'] and params['SKY']['skip']:
        logger.info('Skipping sky background subtraction as requested from parameter file...')
    else:
        logger.info("Starting sky subtraction...")
        try:
            sky.subtract_sky_background(**params['SKY'], in_files=in_files, file_path=params['PATHS']['filePath'])
        except RuntimeError as e:
            raise RuntimeWarning(e)

    # Close reduction
    logger.info("Reduction finished...")
