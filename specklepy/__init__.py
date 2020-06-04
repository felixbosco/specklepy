# Disable automatic IERS download
from astropy.utils import iers
iers.conf.auto_max_age = 60
iers.conf.auto_download = False