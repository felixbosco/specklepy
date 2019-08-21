import os
from configparser import SafeConfigParser
import logging
from logging.config import fileConfig
fileConfig('./config/logging.cfg')


class Reconstructor(object):

    defaults = {}
    config_file_sections = ['Paths', 'Files', 'Parameters']
    parameter_file_keys = ['config_file', 'config_file_name', 'file_name', 'parameter_file', 'parameter_file_name']
    config = './config/reconstructor.cfg'

    def __init__(self, **kwargs):

        # apply default values
        logging.info("Initializing {} with default values...".format(self.__class__.__name__))
        for key in self.defaults:
            setattr(self, key, self.defaults[key])

        # if there is a parameter config file, add it to config
        for key in self.parameter_file_keys:
            if hasattr(self, key):
                self.config = [self.config, getattr(self, key)]

        # read config files
        logging.info("Read config file(s) for {} instance:\n{}".format(self.__class__.__name__, self.config))
        parser = SafeConfigParser()
        parser.optionxform = str  # make option names case sensitive
        found = parser.read(self.config)
        if not found:
            raise ValueError('No config file found!')
        for name in self.config_file_sections:
            self.__dict__.update(parser.items(name))

        # save key word arguments to class attributes
        logging.info("Store key word arguments to {} instance...".format(self.__class__.__name__))
        for key in kwargs:
            setattr(self, key, kwargs[key])

        # evaluating string values of attributes
        logging.info("Transferring str type attributes of {} instance into proper types...".format(self.__class__.__name__))
        for key in self.__dict__:
            try:
                value = getattr(self, key)
                setattr(self, key, eval(value))
            except NameError as e:
                # print(e)
                # value is str type and should be
                pass
            except TypeError as e:
                # print(e)
                # value is already non-str TypeError
                pass
            except SyntaxError as e:
                # print(e)
                # ignoring Paths
                pass

        # creating directories
        for dir in ['input', 'tmp', 'output']:
            if not os.path.isdir(getattr(self, dir)):
                logging.info('Creating {} directory...'.format(dir))
                os.makedirs(getattr(self, dir))


    def __str__(self):
        s = "Class {} object".format(self.__class__.__name__)
        for key in self.__dict__:
            s += "\n{}: {}".format(key, getattr(self, key))
        return s


    def __call__(self, **kwargs):
        logging.info("Starting the reconstruction with a {}...".format(self.__class__.__name__))
        return self.reconstruct(**kwargs)


    def reconstruct(self):
        pass
