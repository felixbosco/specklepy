import os
from logging.config import fileConfig

config_file = os.path.join(os.path.dirname(__file__), 'config/logging.cfg')
fileConfig(config_file)
