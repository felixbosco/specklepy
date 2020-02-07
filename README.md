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
git clone https://github.com/felixbosco/specklepy.git
```

It is recommended to install Specklepy with `pip`:
```bash
pip install .
```

If `pip` is not available to you, you can also execute the setup script. This should work in the same way but has not been tested extensively.
```bash
python setup.py install
```
This installation also creates all the binary scripts that are described below. You may want to double check whether the installation was successful by calling one of the scripts with the `--help` flag:
```bash
holography -h
```

For updating your Specklepy to the latest version, just `git pull` and repeat `pip install .` or `python setup.py install`.

[(top)](table-of-contents)

## Data reduction
The data reduction with Specklepy is divided in two steps: Setup and execution. 
For setting up the files table and parameter file, execute the `setup_reduction` script. 
This gathers files located at `--path` and saves the list to the `--outfile`. 
Note that this function already inspects the fits headers for information on the respective observational setup.
For reading in the file headers properly you have to provide the `--instrument`. 
If this step is not working properly, you can adopt the configuration file `specklepy/config/instruments/cfg`.

The actual reduction is executed by calling the reduction script.
```bash
specklepy_reduction -p your_parameter_file.par
```

[(top)](table-of-contents)


## Image reconstruction
The image reconstruction module is the core of Specklepy, thus referred to as `core`.
Content ...

[(top)](table-of-contents)


## Synthetic image generator
The `synthetic` module is capable of creating synthetic observations from parameter files and reference PSF data, especially short-exposure 'speckle' PSFs. 
The principle code is borrowed from VegaPy, the Virtual Exposure Generator for Astronomy in Python, but the code received multiple accelerations and features. 
A strong focus of the changes with respect to VegaPy was put on readability of the code and documentation. 

For generating exposures just execute the binary:
```bash
generate_exposures -p your_parameter_file.par
```

The parameter files follow `INI` structure syntax and is interpreted by the package parameter file reader. 
Units should be provided by adding `*u.` since the code is translated into `astropy.quantities`, using the `astropy.units` package. 
An example could look like this:
``` 
[TARGET]
band = 'H'
star_table = specklepy/tests/files/example_star_table_arcsec.dat
sky_background = 14.4

[TELESCOPE]
diameter = 8.2*u.m
central_obscuration = 0.14
name = VLT Unit Telescope
psf_source = AiryDisk
psf_resolution = 20.5013*u.mas 
radius = 50.635494509195034*u.mas 


[DETECTOR]
shape = (1024, 1024)
pixel_scale = 0.0106*u.arcsec
dark_current = 0.1*u.electron/u.s
readout_noise = 35*u.electron
system_gain = 17*u.electron/u.adu
optics_transmission = 0.9
quantum_efficiency = 0.9*u.electron/u.ph
# saturation_level = 60000*u.electron

[KWARGS]
DIT = 800*0.2*u.s
nframes = 1
nframes_limit = 1
outfile = specklepy/tests/files/synthetic/airy_200ms.fits
time_stamp = None
dithers = [(0, 0)]
cards={'OBJECT': 'SYNTHETIC', 'OBSTYPE': 'SCIENCE'}
```

[(top)](table-of-contents)


## Development:
These are the current To-Do items:
* General:
  * ~~Make logger write a copy to a .log file~~
  * Create general script with signature "specklepy holography -f ..."
* data reduction module:
  * The `setup_reduction` script should also create a dummy parameter file
  * Implement flatfield correction from flat frames
  * Implement sky subtraction from sky frames
  * Implement linearity correction (wait for D. Thompson)
* core module:
  * Measure "PSF quality" and implement a weighting scheme for frames 
  * Apertures need to throw errors, when the aperture is touching the image edge! This may also be solved by masking the missing data
  * holography - implement secondary source subtraction
  * Tune star finder routines
  * Implement frame save mode, which saves the Fourier transformed images and PSFs to files and then just sums up within a tmp file. This is necessary for "bootstrap resampling", which can deliver uncertainties.
  * **low priority:** implement 'valid' mode, that yields a reconstruction only in the intersection of all file field of views. Current mode is a 'same' mode, which covers the field of the one reference image. We could also think about a 'full' mode that covers the complete observed field, which might be noisy in the dithered outskirts though
  * **low priority:** implement PSF variation by "interpolation" of the reference PSFs to a fraction of the image
* synthetic exposures module:
  * Study whether it makes sense to sample the speckle PSF down before convolution
 
 [(top)](table-of-contents)
 

## Version history
These are the incremental updates between versions:

- *0.4.5dev*: Debugging the scripts and the logging scheme
- *0.4.4*: The setup is now also installing the sub-packages, config files and scripts
- *0.4.3*: The PSF extraction module now contains a mode for creating an effective PSF (in the sense of Anderson & King, 2000), which is oversampling the PSF grid.
- *0.4.2dev*: Developing error propagation of the PSF estimate
- *0.4.1dev*: Developing the data `reduction` module, beginning with sky subtraction
- *0.4.0*: Specklepy received the `synthetic` module for creating synthetic observations
- *0.3.0*: The image reconstruction module now handles cube alignment. This enabled the first successful holographic image reconstruction of multiple dithered data cubes. Note that this bit of code is labeled vesion 0.0.3 in the setup.py.
- *0.2.0*: The holographic image reconstruction is nowfunctional. Note that this bit of code is labeled vesion 0.0.2 in the setup.py.
- *0.1.0*: Note that this bit of code is labeled vesion 0.0.1 in the setup.py.
- *0.0.1dev*

[(top)](table-of-contents)
