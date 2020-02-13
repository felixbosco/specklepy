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

It is recommended to install Specklepy with `pip`. If `pip` is not available to you, you can also execute the setup script. This should work in the same way but has not been tested extensively.
```bash
pip install .
python setup.py install
```

This installation also creates all the binary scripts that are described below. You may want to double check whether the installation was successful by calling one of the scripts with the `--help` flag:
```bash
specklepy_holography -h
```

For updating your Specklepy to the latest version, just `git pull` and repeat `pip install .` or `python setup.py install`.

[(top)](#table-of-contents)



## Data reduction
The data reduction with Specklepy is divided in two steps: Setup and execution.
For setting up the files table and parameter file, execute the `specklepy_reduction_setup` script.
```bash
specklepy_reduction_setup -i your_instrument -p path/to/your/files/ -o files.tab
```

This gathers files located at `--path` and saves the list to the `--outfile`.
Note that this function already inspects the fits headers for information on the respective observational setup.
For reading in the file headers properly you have to provide the `--instrument`.
If this step is not working properly, you can adopt the configuration file `specklepy/config/instruments/cfg`.

The actual reduction is executed by calling the reduction script.
```bash
specklepy_reduction -p your_parameter_file.par
```

[(top)](#table-of-contents)



## Image reconstruction
The image reconstruction module is the core of Specklepy, thus referred to as `core`.
Available reconstructions are the simple shift-and-add (SSA) algorithm and speckle holography.

#### SSA reconstruction
The SSA reconstruction does not take many parameters but reconstructs an image from the input files, specified with the `--file` flag. Temporary files are saved to the `--tmpDir` directory and the reconstruction is saved to the `--outfile`.
```
specklepy_ssa -f your_files_*.fits -t tmp/ -o your_files_ssa.fits
```

#### Holographic reconstruction
For a holographic reconstruction, you can gather all your parameters in a parameter file and execute the script from your working directory:
```
specklepy_holography -p your_parameter_file.par
```

An example parameter file is following `INI` structure. Note that the inDir is guiding to a path and the bash `*` may be used to collect all files matching to the description.
```
[PATHS]
inDir = your_input_files_*.fits # May contain the bash *
outFile = holographic_reconstruction.fits
tmpDir = tmp/
allStarsFile = all_sources.dat
refSourceFile = ref_sources.dat

[STARFINDER]
starfinderFwhm = 5 # pixels
signalToNoiseThreshold = 5  # multiples of sigma

[PSFEXTRACTION]
psfRadius = 45
noiseReferenceMargin = 3
noiseThreshold = 1  # multiples of standard deviation

[APODIZATION]
# apodizationType = Airy # Gaussian or Airy
# apodizationWidth = 4.777 # in pixels
apodizationType = Gaussian
apodizationWidth = 1.645
```

[(top)](#table-of-contents)



## Synthetic image generator
The `synthetic` module is capable of creating synthetic observations from parameter files and reference PSF data, especially short-exposure 'speckle' PSFs.
The principle code is borrowed from VegaPy, the Virtual Exposure Generator for Astronomy in Python, but the code received multiple accelerations and features.
A strong focus of the changes with respect to VegaPy was put on readability of the code and documentation.

For generating exposures just execute the binary:
```bash
generate_exposures -p your_parameter_file.par
```

The parameter files follow `INI` structure syntax and are interpreted by the package parameter file reader.
Units should be provided by adding `*u.` since the code is translated into `astropy.quantities`, using the `astropy.units` package.
If no unit is provided, the objects usually interpret the values in default units.
In the example below, the `sky_background = 14.4` will be interpreted in units mag / arcsec^2.
We refer to the doc-strings in the code for further details.
A parameter file could look like this:
```
[TARGET]
band = 'H'
star_table = your_star_table.dat
sky_background = 14.4

[TELESCOPE]
diameter = 8.2*u.m
central_obscuration = 0.14
name = VLT Unit Telescope
psf_source = AiryDisk
psf_resolution = 20.5*u.mas
radius = 50.6*u.mas

[DETECTOR]
shape = (1024, 1024)
pixel_scale = 0.01*u.arcsec
dark_current = 0.1*u.electron/u.s
readout_noise = 35*u.electron
system_gain = 17*u.electron/u.adu
optics_transmission = 0.9
quantum_efficiency = 0.9*u.electron/u.ph
saturation_level = 60000*u.electron

[KWARGS]
DIT = 1.2*u.s
nframes = 100
nframes_limit = 100
outfile = airy_1200ms.fits
time_stamp = None
dithers = [(0, 0)]
cards={'OBJECT': 'SYNTHETIC', 'OBSTYPE': 'SCIENCE'}
```

[(top)](#table-of-contents)



## Development:
These are the current To-Do items:
* General:
  * **low priority:** Create general script with signature "specklepy holography -f ..."
  * **low priority:** Make command line arguments to positional arguments for required args
* data reduction module:
  * The `specklepy_reduction_setup` script should also create a dummy parameter file
  * Flatfield correction from flat frames
  * Sky subtraction from sky frames
  * Linearity correction (wait for D. Thompson)
* core module:
  * Division of the field of view into chunks with separate PSFs for consideration of variable PSFs.
  * Propagation of uncertainties if provided
  * Measure "PSF quality" and implement a weighting scheme for frames
  * Reconstruction modes that return reconstructions that cover the 'same' field of view as a reference, the 'full' field of view covered by at least one exposure or (**low priority:**) 'valid' field of view that is covered by all exposures
  * Secondary source subtraction for holography function 
  * Bootstrap resampling for estimating the photo- and astrometric uncertainties
  * Add a call of starfinder after the last iteration of the scripts to save the results
  * **low priority:** Apertures need to throw errors, when the aperture is touching the image edge! This may also be solved by masking the missing data
  * **low priority:** Tune star finder routines
  * **low priority:** Implement frame save mode, which saves the Fourier transformed images and PSFs to files and then just sums up within a tmp file. This is necessary for "bootstrap resampling", which can deliver uncertainties.
* synthetic exposures module:
  * Study whether it makes sense to sample the speckle PSF down before convolution

 [(top)](#table-of-contents)



## Version history
These are the incremental updates between versions:

- *0.5.2dev*: Implement noise propagation during data reduction
- *0.5.1*: Implementation of the 'full' mode for reconstructions
- *0.5.0*: The code is now easily installable, creates binary scripts to be used everywhere in the system. Extensive updates on the documentation.
- *0.4.5*: Debugging the scripts and the logging scheme
- *0.4.4*: The setup is now also installing the sub-packages, config files and scripts
- *0.4.3*: The PSF extraction module now contains a mode for creating an effective PSF (in the sense of Anderson & King, 2000), which is oversampling the PSF grid.
- *0.4.2*: Developing error propagation of the PSF estimate
- *0.4.1*: Developing the data `reduction` module, beginning with sky subtraction
- *0.4.0*: Specklepy received the `synthetic` module for creating synthetic observations
- *0.3.0*: The image reconstruction module now handles cube alignment. This enabled the first successful holographic image reconstruction of multiple dithered data cubes. Note that this bit of code is labeled vesion 0.0.3 in the setup.py.
- *0.2.0*: The holographic image reconstruction is nowfunctional. Note that this bit of code is labeled vesion 0.0.2 in the setup.py.
- *0.1.0*: Note that this bit of code is labeled vesion 0.0.1 in the setup.py.
- *0.0.1*: Setup of code structure.

[(top)](#table-of-contents)
