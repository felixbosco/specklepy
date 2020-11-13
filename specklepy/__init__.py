from astropy.utils import iers


# Disable automatic IERS download
iers.conf.auto_max_age = 60
iers.conf.auto_download = False
