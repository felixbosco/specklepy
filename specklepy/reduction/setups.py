import numpy as np
import string
from astropy.table import Table

from specklepy.logging import logging
from specklepy.exceptions import SpecklepyTypeError



def identify_setups(filelist, keywords, return_setups=False):
    """Identify distinct observational setups in a list of files.

    This function identifies distinct observational setups.

    Args:
        filelist (astropy.table.Table):
        keywords (list of str):
        return_setups (bool): Set True to return also the dict of distinct setups. Default is False.
    """

    # Check input parameters
    if not isinstance(filelist, Table):
        raise SpecklepyTypeError('identify_setups', 'filelist', type(filelist), 'astropy.table.Table')

    if not isinstance(keywords, list):
        raise SpecklepyTypeError('identify_setups', 'filelist', type(filelist), 'list')

    if not isinstance(return_setups, bool):
        raise SpecklepyTypeError('identify_setups', 'return_setups', type(return_setups), 'bool')


    # Identifying setups key-by-key
    logging.info("Identifying distinct observational setups in the file list...")
    filelist['Setup'] = [None] * len(filelist['FILE'])

    for key in keywords:
        unique = np.unique(filelist[key].data)
        logging.info("Identified {} setups by keyword {}:".format(len(unique), key))
        logging.info("\t{}".format(unique))
        if len(unique) == 1:
            continue

        for index, setup in enumerate(unique):
            for row in filelist:
                if row[key] == setup:
                    if row['Setup'] is None:
                        row['Setup'] = [str(index)]
                    else:
                        row['Setup'].append(str(index))

    for row in filelist:
        row['Setup'] = ''.join(row['Setup'])

    # Overwrite setup keys by length-1 string
    combinations = np.unique(filelist['Setup'].data)
    for index, combination in enumerate(combinations):
        row_indizes = np.where(filelist['Setup'].data == combination)
        filelist['Setup'][row_indizes] = string.ascii_uppercase[index]


    # Return
    if return_setups:
        # Create setups dict
        setups_dict = {}
        return filelist, setups_dict
    else:
        return filelist