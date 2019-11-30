import os
from configparser import ConfigParser
import astropy.units as u

from specklepy.logging import logging
from specklepy.synthetic.target import Target
from specklepy.synthetic.telescope import Telescope
from specklepy.synthetic.detector import Detector



def get_objects(parameterfile, debug=False):
    """Get objects from parameter file.

    Args:
        parameterfile (str): File from which the objects are instantiated.

    Returns:
        objects (sequence)
    """

    # Check whether files exist
    if not os.path.isfile(parameterfile):
        raise FileNotFoundError("Parameter file {} not found!".format(parameterfile))

    # Prepare objects list
    objects = []

    # Read parameter_file
    parser = ConfigParser(inline_comment_prefixes="#")
    parser.optionxform = str  # make option names case sensitive
    logging.info("Reading parameter file {}".format(parameterfile))
    parser.read(parameterfile)
    for section in parser.sections():
        kwargs = {}
        for key in parser[section]:
            value = parser[section][key]
            try:
                kwargs[key] = eval(value)
            except:
                kwargs[key] = value
            if debug:
                print(key, type(kwargs[key]), kwargs[key])

        if section.lower() == 'target':
            objects.append(Target(**kwargs))
        elif section.lower() == 'telescope':
            objects.append(Telescope(**kwargs))
        elif section.lower() == 'detector':
            objects.append(Detector(**kwargs))

    return objects
