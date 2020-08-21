import configparser
import os
import sys
import yaml

from specklepy.logging import logger


def read(par_file):

    logger.info(f"Reading parameter file {par_file}")

    # Identify type of file
    root, ext = os.path.splitext(par_file)

    if ext == '.yaml':
        return read_yaml(par_file)
    else:
        return read_ini(par_file)


def read_ini(par_file):
    parser = configparser.ConfigParser(inline_comment_prefixes='#')
    parser.optionxform = str  # make option names case sensitive

    # Read in config files
    parser.read(par_file)  # Overwrite defaults

    # Transforming parser information into dict type
    config = {}
    for section in parser.sections():
        config[section] = dict(parser[section])
    return config


def read_yaml(par_file):
    with open(par_file, "r") as yaml_file:
        try:
            config = yaml.load(yaml_file, Loader=yaml.loader.FullLoader)
        except yaml.parser.ParserError as e:
            sys.tracebacklimit = 0
            raise e
    return config


def update_from_file(params, par_file):

    # Read config parameters to update from
    logger.info(f"Updating from parameter file {par_file}")
    update = read(par_file=par_file)

    # Overwrite entries in the input dictionary
    for key in update.keys():
        if isinstance(update[key], dict):
            for kkey in update[key].keys():
                params[key][kkey] = update[key][kkey]
        else:
            params[key] = update[key]

    return params
