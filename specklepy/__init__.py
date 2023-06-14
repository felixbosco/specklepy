from astropy.utils import iers

# from specklepy.plotting.utils import imshow


# Disable automatic IERS download
iers.conf.auto_max_age = 60
iers.conf.auto_download = False
