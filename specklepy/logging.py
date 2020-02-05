import os
import logging.config

config_file = os.path.join(os.path.dirname(__file__), 'config/logging.cfg')
logging.config.fileConfig(config_file)
logging = logging.getLogger('dev')
logging.info("Logging to 'specklepy.log'")
