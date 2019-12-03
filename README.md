# Specklepy - A versatile tool for analyzing astronomical short-exposure ("speckle") data

Specklepy is a versatile tool for data reduction and image reconstruction of short-exposure ("speckle") data. It also contains a module for generating synthetic observations, which is derived from the VegaPy package.


## Table of contents
- [To Do](#to-do)
- [Synthetic image generator](#synthetic-image-generator)
- [Data reduction](#data-reduction)
- [Image reconstruction](#image-reconstruction)
- [Version history](#version-history)


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


## Image reconstruction
The image reconstruction module is the core of Specklepy, thus referred to as `core`.
Content ...


## Data reduction
To be implemented ...


## Synthetic image generator
The `synthetic` module is capable of creating synthetic observations from parameter files and reference PSF data, especially short-exposure 'speckle' PSFs. The principle code is borrowed from VegaPy, the Virtual Exposure Generator for Astronomy in Python, but the code received multiple accelerations and features. A stron focus of the changes with respect to VegaPy was put on readability of the code and documentation.


## Version history
These are the incremental updates between versions:

- *0.4.0*: Specklepy received the `synthetic` module for creating synthetic observations
- *0.3.0*: The image reconstruction module now handles cube alignment. This enabled the first successful holographic image reconstruction of multiple dithered data cubes. Note that this bit of code is labeled vesion 0.0.3 in the setup.py.
- *0.2.0*: The holographic image reconstruction is nowfunctional. Note that this bit of code is labeled vesion 0.0.2 in the setup.py.
- *0.1.0*: Note that this bit of code is labeled vesion 0.0.1 in the setup.py.
- *0.0.1dev*

