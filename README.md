# Specklepy - A versatile tool for analyzing astronomical short-exposure ("speckle") data

Specklepy is a versatile tool for data reduction and image reconstruction of short-exposure ("speckle") data. It also contains a module for generating synthetic observations, which is derived from the VegaPy package.

## Table of contents
[To Do](#to-do)
[Data reduction](#data-reduction)
[Image reconstruction](#image-reconstruction)
[Synthetic image generator](#synthetic-image-generator)


## To Do:
* General:
  * Make logger write a copy to a .log file
  * Implement RAM save mode, which saves the Fourier transformed images and PSFs to files and then just sums up within a tmp file
  * Create general script with signature "specklepy holography -f ..."
* core module:
  * holography - implement secondary source subtraction
  * Tune star finder routines
  * [low priority] implement 'valid' mode, that yields a reconstruction only in the intersection of all file field of views
* data reduction module:
  * Implement flatfield correction
  * Implement linearity correction (wait for D. Thompson)
* synthetic exposures module:
  * Implement dithering
  * Center star coordinates in target around (0, 0), which is then mapped to the center of the detector


## Data reduction
Content ...


## Image reconstruction
Content ...


## Synthetic image generator
Content ...
