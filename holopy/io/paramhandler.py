import os
from configparser import ConfigParser

from holopy.logging import logging
from holopy.utils.listing import Listing


class ParamHandler(object):

    def __init__(self, parameter_file=None, defaults_file=None, essential_attributes=[], make_dirs=[], type_dict={}):
        # Store file names
        self.parameter_file = parameter_file
        self.defaults_file = defaults_file

        # Check whether files exist
        if not os.path.isfile(self.parameter_file):
            raise FileNotFoundError("Parameter file {} not found!".format(self.parameter_file))
        if not os.path.isfile(self.defaults_file):
            raise FileNotFoundError("Defaults file {} not found!".format(self.defaults_file))

        # Read parameter_file
        parser = ConfigParser()
        parser.optionxform = str  # make option names case sensitive
        logging.info("Reading parameter file {}".format(self.parameter_file))
        parser.read(self.parameter_file)
        for section in parser.sections():
            for key in parser[section]:
                value = parser[section][key]
                # Interprete data type
                try:
                    setattr(self, key, eval(value))
                except:
                    setattr(self, key, value)


        # Complete attribute list from defaults file
        defaults = ConfigParser()
        defaults.optionxform = str  # make option names case sensitive
        logging.info("Reading defaults file {}".format(self.defaults_file))
        defaults.read(self.defaults_file)
        for attr in essential_attributes:
            if not hasattr(self, attr):
                attr_set = False
                for section in defaults.sections():
                    for key in defaults[section]:
                        if key == attr:
                            value = defaults[section][key]
                            # Interprete data type
                            try:
                                setattr(self, key, eval(value))
                            except:
                                setattr(self, key, value)
                            attr_set = True
                if not attr_set:
                    logging.warning("Essential parameter <{}> not found in parameter file or config file!".format(attr))

        self.makedirs(dir_list=make_dirs)


    def makedirs(self, dir_list):
        """
        This function makes sure that the paths exist and creates if not.
        """
        for key in dir_list:
            path = getattr(self, key)
            if not os.path.exists(path):
                logging.info('Creating {} directory {}'.format(key, path))
                os.makedirs(path)
