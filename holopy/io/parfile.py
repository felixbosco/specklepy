from configparser import SafeConfigParser

from holopy.logging import logging


class Parfile(object):

    def __init__(self, filename=None):
        if filename is None:
            return 0

        else:
            parser = SafeConfigParser()
            parser.optionxform = str  # make option names case sensitive
            found = parser.read(filename)
