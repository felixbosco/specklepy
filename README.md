# Specklepy - A versatile tool for analyzing astronomical short-exposure ("speckle") data

Specklepy is a versatile tool for data reduction and image reconstruction of short-exposure ("speckle") data. It also contains a module for generating synthetic observations, which is derived from the VegaPy package.


## Table of contents
- [Installation](#installation)
- [Data reduction](#data-reduction)
- [Image reconstruction](#image-reconstruction)
- [Synthetic image generator](#synthetic-image-generator)
- [Development](#development)
- [Version history](#version-history)

## Installation
The source code of Specklepy is available at github and can be downloaded by cloning the git repository to an arbitrary directory. 
`cd` to the directory of choice and execute:
```bash
git clone 
```

It is recommended to install Specklepy with `pip`:
```bash
pip install .
```

If `pip` is not available to you, you can also execute the setup script:
```bash
python setup.py install
```

This installation also creates all the binary scripts that are described below. You may want to double check whether the installation was successful by calling one of the scripts with the `--help` flag:
 ```bash
holography -h
```

## Data reduction
To be implemented ...


## Image reconstruction
The image reconstruction module is the core of Specklepy, thus referred to as `core`.
Content ...


## Synthetic image generator
The `synthetic` module is capable of creating synthetic observations from parameter files and reference PSF data, especially short-exposure 'speckle' PSFs. The principle code is borrowed from VegaPy, the Virtual Exposure Generator for Astronomy in Python, but the code received multiple accelerations and features. A stron focus of the changes with respect to VegaPy was put on readability of the code and documentation.


## Development:
These are the current To-Do items:
* General:
  * Make logger write a copy to a .log file
  * Create general script with signature "specklepy holography -f ..."
* data reduction module:
  * Implement flatfield correction from flat frames
  * Implement sky subtraction from sky frames
  * Implement linearity correction (wait for D. Thompson)
* core module:
  * Measure "PSF quality" and implement a weighting scheme for frames 
  * Apertures need to throw errors, when the aperture is touching the image edge! This may also be solved by masking the missing data
  * holography - implement secondary source subtraction
  * Tune star finder routines
  * Implement frame save mode, which saves the Fourier transformed images and PSFs to files and then just sums up within a tmp file. This is necessary for "bootstrap resampling", which can deliver uncertainties.
  * [low priority] implement 'valid' mode, that yields a reconstruction only in the intersection of all file field of views. Current mode is a 'same' mode, which covers the field of the one reference image. We could also think about a 'full' mode that covers the complete observed field, which might be noisy in the dithered outskirts though
  * [low priority] implement PSF variation by "interpolation" of the reference PSFs to a fraction of the image
* synthetic exposures module:
  * Study whether it makes sense to sample the speckle PSF down before convolution
  

## Version history
These are the incremental updates between versions:

- *0.4.4dev*: Debugging the scripts and the logging scheme
- *0.4.4*: The setup is now also installing the sub-packages, config files and scripts
- *0.4.3*: The PSF extraction module now contains a mode for creating an effective PSF (in the sense of Anderson & King, 2000), which is oversampling the PSF grid.
- *0.4.2dev*: Developing error propagation of the PSF estimate
- *0.4.1dev*: Developing the data `reduction` module, beginning with sky subtraction
- *0.4.0*: Specklepy received the `synthetic` module for creating synthetic observations
- *0.3.0*: The image reconstruction module now handles cube alignment. This enabled the first successful holographic image reconstruction of multiple dithered data cubes. Note that this bit of code is labeled vesion 0.0.3 in the setup.py.
- *0.2.0*: The holographic image reconstruction is nowfunctional. Note that this bit of code is labeled vesion 0.0.2 in the setup.py.
- *0.1.0*: Note that this bit of code is labeled vesion 0.0.1 in the setup.py.
- *0.0.1dev*
