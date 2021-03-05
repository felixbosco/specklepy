import glob
from IPython import embed
import os

from astropy.io import fits
from astropy.table import Table

from specklepy.io import config
from specklepy.io.filearchive import ReductionFileArchive
from specklepy.logging import logger
from specklepy.reduction import flat, sky


def inspect(files, keywords, save=None, recursive=False, debug=False):

    # Find files
    if len(files) == 1:
        path = files[0]
        if '*' in path:
            files = glob.glob(path, recursive=recursive)
        else:
            files = glob.glob(os.path.join(path, '*fits'), recursive=recursive)

    # Terminal output
    if len(files):
        logger.info(f"Inspecting {len(files)} files for the FITS header keywords: {keywords}")
        files.sort()
    else:
        logger.error(f"Found no files in {files}!")
        raise RuntimeError(f"Found no files in {files}!")

    # Initialize output table
    out_table = Table(names=['FILE'] + keywords, dtype=[str] * (1 + len(keywords)))

    # Iterate through files and store header keywords
    for file in files:
        hdr = fits.getheader(filename=file)
        row = [file]
        for keyword in keywords:
            try:
                row.append(str(hdr[keyword]))
            except KeyError:
                row.append('--')
        out_table.add_row(row)

    # Store the table if requested
    if save is not None:
        out_table.write(save, format='ascii.fixed_width')

    # Present table
    print(out_table)


def setup(path, instrument, par_file=None, list_file=None, sort_by=None, recursive=False):
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
        recursive (bool, optional):
            Search for files in a recursive way, that is all sub-directories.
    """

    # Defaults
    default_cards = ['OBSTYPE', 'OBJECT', 'FILTER', 'EXPTIME', 'DIT', 'nFRAMES', 'DATE', 'SUBWIN']
    dtypes = [str, str, str, float, float, int, str, str]
    instrument_config_file = os.path.join(os.path.dirname(__file__), '../config/instruments.cfg')

    # Read config
    configs = config.read(instrument_config_file)
    instrument_cards = configs[instrument.upper()]

    # Double check whether all aliases are defined
    cards = []
    for card in default_cards:
        try:
            cards.append(instrument_cards[card])
        except KeyError:
            logger.info(
                f"Dropping header card {card} from setup identification, as there is no description in the config file."
                f"\nCheck out {instrument_config_file} for details.")
            cards.append(None)
    for card, dtype, header_card in zip(cards, dtypes, default_cards):
        if card is None:
            cards.remove(card)
            dtypes.remove(dtype)
            default_cards.remove(header_card)

    # Apply fall back values
    if path is None:
        path = '.'
    if list_file is None:
        list_file = 'files.tab'
    if par_file is None:
        par_file = 'reduction.yaml'

    # Find files
    if '*' in path:
        files = glob.glob(path, recursive=recursive)
    else:
        files = glob.glob(os.path.join(path, '*fits'), recursive=recursive)
    if len(files):
        logger.info(f"Found {len(files)} file(s)")
        files.sort()
    else:
        logger.error(f"Found no files in {path}!")
        raise RuntimeError(f"Found no files in {path}!")

    # Initialize a file archive
    raw_files = ReductionFileArchive(files, cards=cards, dtypes=dtypes, names=default_cards, sort_by=sort_by)
    raw_files.identify_setups(['FILTER', 'EXPTIME'])
    raw_files.write_table(file_name=list_file)

    # Write dummy parameter file for the reduction
    _, ext = os.path.splitext(par_file)
    if 'yaml' in ext:
        logger.info(f"Creating default reduction YAML parameter file {par_file}")
        par_file_content = f"PATHS:\n  filePath: {raw_files.in_dir}\n  fileList: {list_file}\n  outDir: Science/" \
                           f"\n  tmpDir: Master//\n  prefix: r" \
                           f"\n\nDARK:\n  masterDarkFile: MasterDark.fits" \
                           f"\n\nFLAT:\n  masterFlatFile: MasterFlat.fits" \
                           f"\n\nSKY:\n  method: scalar"
    else:
        logger.info(f"Creating default reduction INI parameter file {par_file}")
        par_file_content = f"[PATHS]\nfilePath = {raw_files.in_dir}\nfileList = {list_file}\noutDir = Science/" \
                           f"\ntmpDir = Master/\nprefix = r" \
                           f"\n\n[DARK]\nmasterDarkFile = MasterDark.fits" \
                           f"\n\n[FLAT]\nmasterFlatFile = MasterFlat.fits" \
                           f"\n\n[SKY]\nmethod = scalar"
    with open(par_file, 'w+') as par_file:
        par_file.write(par_file_content)


def full_reduction(params, debug=False):
    """Execute a full reduction following the parameters in the `params` dictionary.

    TODO: Split this function into the parts and sort into the other modules

    Args:
        params (dict):
            Dictionary with all the settings for reduction.
        debug (bool, optional):
            Show debugging information.
    """

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')

    # Extract dictionaries from params
    paths = params.get('PATHS')
    flat_fielding = params.get('FLAT')
    sky_subtraction = params.get('SKY')
    options = params.get('OPTIONS')


    # (0) Read file list table
    logger.info("Reading file list ...")
    in_files = ReductionFileArchive(file_list=paths.get('fileList'), in_dir=paths.get('filePath'),
                                    out_dir=paths.get('outDir'))
    logger.info('\n' + str(in_files.table))

    # (1) Initialize directories and reduction files
    if not os.path.isdir(paths.get('outDir')):
        logger.debug(f"Making directory {paths.get('outDir')}")
        os.makedirs(paths.get('outDir'))
    if not os.path.isdir(paths.get('tmpDir')):
        logger.debug(f"Making directory {paths.get('tmpDir')}")
        os.makedirs(paths.get('tmpDir'))
    product_files, is_product_file = in_files.make_product_file_names(prefix=paths.get('prefix'),
                                                                      return_table_mask=True)

    # (2) Flat fielding
    if flat_fielding.get('skip', False):
        logger.info('Skipping flat fielding as requested from parameter file...')
    else:
        flat_files = in_files.get_flats()
        if len(flat_files) == 0:
            logger.warning("Did not find any flat field observations. No flat field correction will be applied!")
        else:
            logger.info("Starting flat field correction...")
            master_flat = flat.MasterFlat(flat_files, file_name=flat_fielding.get('masterFlatFile'),
                                          file_path=paths.get('filePath'), out_dir=paths.get('tmpDir'),
                                          new=not flat_fielding.get('reuse', False))
            if not flat_fielding.get('reuse', False):
                master_flat.combine(method=flat_fielding.get('method'))

            for index in range(len(product_files)):
                product_file = in_files.initialize_product_file(index=index)
                sub_window = in_files.table['SUBWIN'][is_product_file][index]
                master_flat.run_correction(file_list=product_file, file_path=None, sub_windows=sub_window,
                                           full_window=options.get('full_window'))

    # (3) Linearization
    # TODO: Implement linearization

    # (4) Sky subtraction
    if sky_subtraction.get('skip', False):
        logger.info('Skipping sky background subtraction as requested from parameter file...')
    else:
        logger.info("Starting sky subtraction...")
        sky.subtract_sky_background(**sky_subtraction, in_files=in_files, out_files=product_files,
                                    file_path=paths.get('filePath'), tmp_dir=paths.get('tmpDir'))

    # Close reduction
    logger.info("Reduction finished!")
