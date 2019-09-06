import os
from configparser import ConfigParser

from holopy.logging import logging
from holopy.utils.listing import Listing


class ParamHandler(object):

    def __init__(self, parameter_file=None, defaults_file=None, essential_attributes=[]):
        # Store file names
        self.parameter_file = parameter_file
        self.defaults_file = defaults_file

        # Check whether files exist
        if not os.path.isfile(self.parameter_file):
            raise FileNotFoundError
        if not os.path.isfile(self.defaults_file):
            raise FileNotFoundError


        # Read parameter_file
        parser = ConfigParser()
        parser.optionxform = str  # make option names case sensitive
        logging.info("Reading parameter file {}".format(self.parameter_file))
        parser.read(self.parameter_file)
        for section in parser.sections():
            setattr(self, section, Listing())
            for key in parser[section]:
                setattr(getattr(self, section), key, parser[section][key])

        return None
        
        # Complete attribute list from defaults file
        for attr in essential_attributes:
            pass


    def mkdirs(self):
        # creating directories
        # for dir in ['input', 'tmp', 'output']:
        #     if not os.path.isdir(getattr(self, dir)):
        #         logging.info('Creating {} directory...'.format(dir))
        #         os.makedirs(getattr(self, dir))
        pass
