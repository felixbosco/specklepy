import os
from configparser import ConfigParser

from holopy.logging import logging


class ParamHandler(object):

    def __init__(self, parameter_file=None, defaults_file=None, essential_attributes=[]):
        # Store file names
        self.parameter_file = parameter_file
        self.defaults_file = defaults_file

        # Read parameter_file
        parser = ConfigParser()
        parser.optionxform = str  # make option names case sensitive
        logging.info("Reading parameter file {}".format(self.parameter_file))
        with open(self.parameter_file, 'r') as configfile:
            configs = parser.read(configfile)
        for section in configs:#.sections():
            print(section)

        return None
        found = parser.read(parameter_file)
        if not found:
            raise ValueError('Parameter file {} not found!'.format(parameter_file))

        for name in self.config_file_sections:
            self.__dict__.update(parser.items(name))

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
