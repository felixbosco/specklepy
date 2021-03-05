import glob
import os

from specklepy.io import config
from specklepy.io.config import read
from specklepy.io.filearchive import ReductionFileArchive
from specklepy.logging import logger
from specklepy.reduction import dark, flat


class DataReduction(object):

    config_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '../config/reduction.cfg'))

    def __init__(self, paths=None, dark=None, linear=None, flat=None, sky=None, options=None, **kwargs):
        self.paths = paths
        self.dark = dark
        self.linear = linear
        self.flat = flat
        self.sky = sky
        self.options = options
        for attr in ['paths', 'dark', 'linear', 'flat', 'sky', 'options']:
            if getattr(self, attr) is None:
                logger.debug(f"Transforming parameters from section {attr.upper()!r} in the config file to {attr!r}...")
                setattr(self, attr, kwargs.get(attr.upper()))

        # Initialize file archive
        self.files = ReductionFileArchive(file_list=self.paths.get('fileList'), in_dir=self.paths.get('filePath'),
                                          out_dir=self.paths.get('outDir'))

    @classmethod
    def from_file(cls, file_name):
        logger.info(f"Configuring data reduction from config file {cls.config_file!r}")
        config = read(par_file=cls.config_file)
        logger.info(f"Updating data reduction configuration from parameter file {file_name!r}")
        params = read(par_file=file_name)
        for section in params:
            config[section].update(**params[section])
        return cls(**config)

    def run(self):
        self.initialize_directories()
        self.initialize_product_files()
        if not self.dark.get('skip'):
            self.run_dark_correction()
        if not self.linear.get('skip'):
            self.run_linearization()
        if not self.flat.get('skip'):
            self.run_flat_fielding()
        if not self.sky.get('skip'):
            self.run_sky_subtraction()

        # Close reduction
        logger.info("Reduction finished!")

    def initialize_directories(self):
        if not os.path.isdir(self.paths.get('outDir')):
            logger.debug(f"Making directory {self.paths.get('outDir')}")
            os.makedirs(self.paths.get('outDir'))
        if not os.path.isdir(self.paths.get('tmpDir')):
            logger.debug(f"Making directory {self.paths.get('tmpDir')}")
            os.makedirs(self.paths.get('tmpDir'))

    @classmethod
    def set_up(cls, path, instrument, par_file=None, list_file=None, sort_by=None, recursive=False):
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
        raw_files.add_dark_column()
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

    def initialize_product_files(self):
        self.files.add_product_file_column()
        if self.options.get('clearEarlierProductFiles', False):
            logger.info("Removing data product files from earlier reductions...")
            for product_file_path in self.files.product_file_paths:
                try:
                    os.remove(product_file_path)
                except FileNotFoundError:
                    pass

    def run_dark_correction(self):

        # Identify dark setups
        dark_setups = self.files.get_dark_setups()
        master_darks = {}

        # Create a master dark for each setup
        if not self.dark.get('reuse'):
            for setup in dark_setups:
                darks = self.files.filter({'OBSTYPE': 'DARK', 'SETUP': setup})
                master_dark = dark.MasterDark(file_list=darks, file_path=self.files.in_dir, setup=setup,
                                              file_name=self.dark.get('masterDarkFile'),
                                              out_dir=self.paths.get('tmpDir'))
                master_dark.combine(number_frames=self.dark.get('numberFrames'))
                master_dark.write()
                master_darks[setup] = master_dark.path

        # Apply dark subtraction
        for setup in dark_setups:
            master_dark = dark.MasterDark.from_file(master_darks[setup])
            for p, product_file_path in enumerate(self.files.product_file_paths):
                # Skip files not allocated
                if self.files.table['DARK'][p] != setup:
                    continue

                # Initialize file if not existing
                if not os.path.isfile(product_file_path):
                    self.files.initialize_product_file(index=p)

                # Extract sub-window for file
                sub_window = self.files.table['SUBWIN'].data[p]

                # Subtract master dark from file
                master_dark.subtract(product_file_path, sub_window=sub_window)

    def run_flat_fielding(self):

        # Identify flat filters
        flat_filters = self.files.get_flat_filters()
        master_flats = {}

        # Create a master flat for each filter
        if not self.flat.get('reuse'):
            for filter in flat_filters:
                flats = self.files.filter({'OBSTYPE': 'FLAT', 'FILTER': filter}, namekey='PRODUCT')
                master_flat = flat.MasterFlat(file_list=flats, file_name=self.flat.get('masterFlatFile'),
                                              file_path=self.files.in_dir, out_dir=self.paths.get('tmpDir'),
                                              filter=filter)
                master_flat.combine(method=self.flat.get('method'))
                master_flat.write()
                master_flats[filter] = master_flat.path

        # Apply flat field correction
        for filter in flat_filters:
            master_flat = flat.MasterFlat.from_file(master_flats[filter])
            for p, product_file_path in enumerate(self.files.product_file_paths):
                # Skip files from different filter
                if self.files.table['FILTER'][p] != filter:
                    continue

                # Initialize file if not existing
                if not os.path.isfile(product_file_path):
                    self.files.initialize_product_file(index=p)

                # Extract sub-window for file
                sub_window = self.files.table['SUBWIN'].data[p]

                # Normalize product file with master flat
                master_flat.run_correction(file_list=[product_file_path], sub_windows=[sub_window])

    def run_linearization(self):
        raise NotImplementedError("Linearization is not implemented yet!")

    def run_sky_subtraction(self):
        raise NotImplementedError("Sky subtraction is not implemented yet!")
