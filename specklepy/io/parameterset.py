import os
import configparser

from specklepy.exceptions import SpecklepyTypeError
from specklepy.logging import logger
from specklepy.io.filearchive import FileArchive
from specklepy.mock.target import Target
from specklepy.mock.telescope import Telescope
from specklepy.mock.detector import Detector


class ParameterSet(object):

    """Namespace class for storing parameters.

    This class stores all the parameters. The values are provided by reading in a default/ config file first
    (if a file is provided) and then reading in the actual parameter file to overwrite the default values. The
    parameters are stored hierarchically following the sections of the defaults or parameter files to allow the same
    names in multiple sections. The parameters are stored within Section class instances, see below, and are thereby
    accessible via attribute calls.

    In a second step this class is checking a list of essential parameters for existence and creates directories, if
    the arguments are not None.

    Attributes:
        ...
    """

    def __init__(self, parameter_file, defaults_file=None, essential_attributes=None, make_dirs=None,
                 store_mode='attr'):
        """Create a ParameterSet instance.

        Args:
            parameter_file (str):
                Path to a file that contains parameters
            defaults_file (str, optional):
                Path to a file that contains default parameters.
            essential_attributes (dict, optional):
                Dict of attributes that are essential and are thus filled from
                the defaults file, if not provided in the parameter file.
                Absent parameters cause exceptions.
            make_dirs (list, optional):
                List of directory paths to create, if they are not existing yet.
            store_mode (str, optional):
                Mode of storing the parsed parameters. In 'attr' mode, all values are stored inside Section instances.
                In 'dict' mode, the parameters are stored as nested dictionaries.
        """

        # Check input parameters
        if isinstance(parameter_file, str):
            self.parameter_file = os.path.abspath(parameter_file)
        else:
            raise SpecklepyTypeError('ParameterSet', argname='parameter_file',
                                     argtype=type(parameter_file), expected='str')

        if isinstance(defaults_file, str):
            self.defaults_file = os.path.abspath(defaults_file)
        elif defaults_file is None:
            self.defaults_file = None
        else:
            raise SpecklepyTypeError('ParameterSet', argname='defaults_file',
                                     argtype=type(defaults_file), expected='str')

        if essential_attributes is None:
            essential_attributes = {}
        if make_dirs is None:
            make_dirs = []

        # Set up config parser
        parser = configparser.ConfigParser(inline_comment_prefixes='#')
        parser.optionxform = str  # make option names case sensitive

        # Read in config files
        if self.defaults_file:
            logger.info(f"Reading defaults file {self.defaults_file}")
            parser.read(self.defaults_file)
        logger.info(f"Reading parameter file {self.parameter_file}")
        parser.read(parameter_file)  # Overwrite defaults

        # Store parser attributes
        if store_mode == 'attr':
            # ... into ParameterSet Sections
            for section in parser.sections():
                setattr(self, section.lower(), Section(parser[section]))
        elif store_mode == 'dict':
            # TODO: Implement this properly
            for section in parser.sections():
                setattr(self, section.lower(), dict(parser[section]))
        else:
            raise ValueError(f"Storing mode {store_mode} is unknown!")

        # Check for the essential parameters
        for section in essential_attributes.keys():
            if not hasattr(self, section):
                raise AttributeError(f"ParameterSet does not have essential section {section}!")
            for attr in essential_attributes[section]:
                if not hasattr(getattr(self, section), attr):
                    raise AttributeError(f"ParameterSet.{section} does not have essential attribute {attr}!")

        # Create directories
        self.makedirs(dir_list=make_dirs)

        # Create file lists
        if hasattr(self, 'paths'):
            if hasattr(self.paths, 'inDir'):
                self.inFiles = FileArchive(self.paths.inDir).files
        # try:
        #     self.inFiles = FileManager(self.paths.inDir).files
        # except AttributeError:
        #     logger.warning("ParameterSet instance is not storing 'inFiles' due to missing entry 'inDir' parameter in "
        #                    "parameter file!")

    def makedirs(self, dir_list):
        """
        This function makes sure that the paths exist and creates if not.
        """

        for key in dir_list:
            path = getattr(self.paths, key)
            path = os.path.dirname(path) + '/'  # Cosmetics to allow for generic input for inDir
            if not os.path.exists(path):
                logger.info(f"Creating {key} directory {path}")
                os.makedirs(path)


class Section(object):

    """Namespace class that stores section options.

    Attributes:
        ...
    """

    def __init__(self, options=None):
        """Create a Section instance from an options dict.

        Args:
            options (configparser.SectionProxy, optional):
                All section options to be stored.
        """

        if not isinstance(options, (dict, configparser.SectionProxy)):
            raise SpecklepyTypeError('Section', argname='options', argtype=type(options), expected='dict')

        for key in options.keys():
            value = options[key]
            try:
                setattr(self, key, eval(value))
            except:
                setattr(self, key, value)


class ReductionParameterSet(ParameterSet):

    """ParameterSet sub-class with default parameters.

    This ParameterSet child class is created with some default values for data reduction.

    Attributes:
        See ParameterSet.
    """

    def __init__(self, parfile):
        """Create a ReductionParameterSet instance.

        Args:
            parfile (str):
                Name of the parameter file to create the ParameterSet from.
        """
        defaults_file = os.path.join(os.path.dirname(__file__), '../config/reduction.cfg')
        defaults_file = os.path.abspath(defaults_file)
        essential_attributes = {'paths': ['filePath', 'fileList', 'outDir', 'tmpDir', 'filePrefix'],
                                'setup': ['setupKeywords'],
                                'flat': ['skip', 'masterFlatFile'],
                                'sky': ['skip', 'source', 'method', 'sigmaClip']}

        super().__init__(parameter_file=parfile,
                         defaults_file=defaults_file,
                         essential_attributes=essential_attributes,
                         make_dirs=['tmpDir'])


class HolographyParameterSet(ParameterSet):

    """ParameterSet sub-class with default parameters.

    This ParameterSet child class is created with some default values for the holography algorithm.

    Attributes:
        See ParameterSet.
    """

    def __init__(self, parfile):
        """Create a HolographyParameterSet instance.

        Args:
            parfile (str):
                Name of the parameter file to create the ParameterSet from.
        """

        # Default values
        defaults_file = os.path.join(os.path.dirname(__file__), '../config/holography.cfg')
        defaults_file = os.path.abspath(defaults_file)
        essential_attributes = {'paths': ['inDir', 'tmpDir', 'outFile', 'alignmentReferenceFile', 'refSourceFile'],
                                'starfinder': ['starfinderFwhm', 'noiseThreshold'],
                                'psfextraction': ['mode', 'psfRadius', 'noiseThreshold', 'noiseReferenceMargin',
                                                  'fieldSegmentation'],
                                'apodization': ['apodizationWidth', 'apodizationType'],
                                'options': ['reconstructionMode', 'varianceExtensionName']}

        # Read parameters from file
        super().__init__(parameter_file=parfile,
                         defaults_file=defaults_file,
                         essential_attributes=essential_attributes,
                         make_dirs=['tmpDir'])


class GeneratorParameterSet(ParameterSet):

    def __init__(self, parfile):

        essential_attributes = {'target': ['band'],
                                'telescope': ['diameter', 'psf_source'],
                                'detector': ['shape', 'pixel_scale'],
                                'parameters': ['exposure_time', 'n_frames']}

        super().__init__(parfile, defaults_file=None, essential_attributes=essential_attributes)

        from IPython import embed
        embed()

        self.target = Target(**self.target)
        self.telescope = Telescope(**self.telescope)
        self.detector = Detector(**self.detector)
