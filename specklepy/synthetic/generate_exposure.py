import os
from datetime import datetime
import numpy as np
from configparser import ConfigParser

from astropy.io import fits
from astropy.units import Quantity

from specklepy.logging import logger
from specklepy.synthetic.target import Target
from specklepy.synthetic.telescope import Telescope
from specklepy.synthetic.detector import Detector


def generate_exposure(target, telescope, detector, DIT,
                      nframes=1, nframes_limit=100, dithers=None,
                      outfile='exposure.fits', time_stamp=None,
                      debug=False, **kwargs):
    """Generate synthetic exposures from target, telescope and detector objects.

    The function generate_exposure() is the central function of the synthetic
    module. It takes one instance each  of the classes Target, Telescope, and
    Detector and the discrete exposure time DIT. Then, it creates a number of
    files depending on the number of requested frames ('nframes') and frame
    limit per file ('nframes_limit').
    To distribute the synthetic exposures to multiple files, for instance if the
    size of the individual file would become too large, just set
    'nframes_limit' to a smaller value (default is 100).

    Args:
        target (specklepy.synthetic.Target):
            Target instance that will be 'observed'.
        telescope (specklepy.synthetic.Telescope):
            Telescope instance that will 'observe'.
        detector (specklepy.synthetic.Detector):
            Detector instance that will be 'exposed'.
        DIT (astropy.units.Quantity):
            Discrete integration time for each exposure.
        nframes (int, optional):
            Number of frames that will be generated. Default is 1.
        nframes_limit (int, optional):
            Maximum number of frames per outfile. Default is 100.
        dithers (list, optional):
            List of dither positions, offsets from the Target instance center. Default is None.
        outfile (str, optional):
            Base name of the output file. If multiple files are requested, the index will be appended. Default is
            'exposure.fits'.
        time_stamp (str, optional):
            Time stamp that can be added to outfile. Can be either 'start', 'end' or None to suppress adding the time
            stamp. Default is None.
        debug (bool, optional):
            Show debugging information. Default is False.
    """

    # Compute number of files
    nfiles = nframes // nframes_limit
    nframes_left = nframes % nframes_limit
    if nframes_left:
        logger.info(f"Creating {nfiles} files with {nframes_limit} synthetic exposures and adding {nframes_left} "
                    f"frames to an additional file")
    else:
        logger.info(f"Creating {nfiles} files with {nframes_limit} synthetic exposures")

    # Add a time stamp to file name, if not None
    outdir, outfile = os.path.split(outfile)
    now = datetime.now()
    time_str = now.strftime('%Y%m%d_%H%M%S')
    if time_stamp == 'end':
        outfile = outfile.replace('.fits', f"_{time_str}.fits")
    elif time_stamp == 'start':
        outfile = f"{time_str}_{outfile}"
    elif time_stamp is None:
        pass
    outfile = os.path.join(outdir, outfile)

    # Initialize fits header
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((nframes_limit,) + detector.shape)
    hdu.header.set('DIT', DIT.value, DIT.unit)
    hdu.header.set('DATE', str(datetime.now()))
    if 'cards' in kwargs:
        for key in kwargs['cards']:
            hdu.header.set(key, kwargs['cards'][key])
    # Add object attributes to header information
    skip_attributes = {'target': ['data', 'stars'],
                       'telescope': ['psf', 'psf_frame'],
                       'detector': ['shape', 'array']}
    for object in [target, telescope, detector]:
        dict = object.__dict__
        object_name = object.__name__

        for key in dict:
            if key in skip_attributes[object_name]:
                continue
            card = f"HIERARCH SPECKLEPY {object_name.upper()} {key.upper()}"
            if isinstance(dict[key], Quantity):
                hdu.header.set(card, f"{dict[key].value:.3e}", dict[key].unit)  # Appending the unit of a Quantity to the comment
            elif isinstance(dict[key], str):
                hdu.header.set(card, os.path.basename(dict[key]))  # Suppress (long) relative paths
            elif isinstance(dict[key], tuple):
                _tuple = tuple([x.value for x in dict[key]])  # Separating tuple unit from values
                hdu.header.set(card, str(_tuple), dict[key][0].unit)
            else:
                hdu.header.set(card, dict[key])

    # Write header to one or more files, depending on 'nframes' and 'nframes_limit'
    outfiles = []
    for n in range(nfiles):
        filename = outfile.replace('.fits', f"_{n + 1}.fits")
        outfiles.append(filename)
        hdu.writeto(filename, overwrite=True)
    # Create file for the left over frames
    if nframes_left > 0:
        filename = outfile.replace('.fits', f"_{nfiles + 1}.fits")
        outfiles.append(filename)
        hdu.data = np.zeros((nframes_left,) + detector.shape)
        hdu.writeto(filename, overwrite=True)

    # Initialize parameters for frame computation
    if 'readout_time' in kwargs:
        skip_frames = int(kwargs['readout_time'] / telescope.psf_timestep)
    else:
        skip_frames = 0

    # Computation of frames
    frame_counter = 0
    for outfile_index, outfile in enumerate(outfiles):
        with fits.open(outfile, mode='update') as hdulist:
            # Get a new field of view for each file to enable dithering between files
            if dithers is not None:
                try:
                    dither = dithers[outfile_index]
                except IndexError:
                    raise RuntimeError(f"Expected {len(outfiles)} dither positions but received only {len(dithers)}!")
            else:
                dither = None
            photon_rate_density = target.get_photon_rate_density(FoV=detector.FoV, resolution=telescope.psf_resolution,
                                                                 dither=dither)

            for index in range(hdulist[0].header['NAXIS3']):
                photon_rate = telescope.get_photon_rate(photon_rate_density, integration_time=DIT, debug=debug)
                counts = detector.get_counts(photon_rate=photon_rate,
                                             integration_time=DIT,
                                             photon_rate_resolution=target.resolution, debug=debug)
                counts = counts.decompose()
                hdulist[0].data[index] = counts.value

                try:
                    telescope.psf_frame += skip_frames
                except TypeError:
                    pass

                frame_counter += 1
                print("\rSaving exposure {:4}/{:4} to file {}...".format(frame_counter, nframes, outfile), end='')
    print("")


def get_objects(parameterfile, debug=False):
    """Get objects from parameter file.

    Args:
        parameterfile (str):
            File from which the objects are instantiated.
        debug (bool, optional):
            Show debugging information.

    Returns:
        objects (dict):
            Dict containing the parameter file section as a key word and the objects (or kwargs dict) as values.
    """

    # Check whether files exist
    if not os.path.isfile(parameterfile):
        raise FileNotFoundError(f"Parameter file {parameterfile} not found!")

    # Prepare objects list
    objects = {}

    # Read parameter_file
    parser = ConfigParser(inline_comment_prefixes="#")
    parser.optionxform = str  # make option names case sensitive
    logger.info(f"Reading parameter file {parameterfile}")
    parser.read(parameterfile)
    for section in parser.sections():
        kwargs = {}
        for key in parser[section]:
            value = parser[section][key]
            try:
                kwargs[key] = eval(value)
            except SyntaxError:
                try:
                    kwargs[key] = Quantity(value)
                except:
                    kwargs[key] = value
            except:
                kwargs[key] = value
            if debug:
                print(key, type(kwargs[key]), kwargs[key])

        if section.lower() == 'target':
            objects['target'] = Target(**kwargs)
        elif section.lower() == 'telescope':
            objects['telescope'] = Telescope(**kwargs)
        elif section.lower() == 'detector':
            objects['detector'] = Detector(**kwargs)
        elif section.lower() in ['params', 'parameters', 'kwargs']:
            objects['parameters'] = kwargs

    return objects
