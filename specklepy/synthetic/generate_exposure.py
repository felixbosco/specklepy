import os
from datetime import datetime
import numpy as np
from configparser import ConfigParser
from astropy.io import fits
import astropy.units as u

from specklepy.logging import logger
from specklepy.synthetic.target import Target
from specklepy.synthetic.telescope import Telescope
from specklepy.synthetic.detector import Detector



def generate_exposure(target, telescope, detector, DIT, nframes=1, nframes_limit=100, dithers=None, outdir=None, outfile='exposure.fits', time_stamp='end', debug=False, **kwargs):
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
        telescope (specklepy.synthetic.Telescope):
        detector (specklepy.synthetic.Detector):
        DIT (astropy.units.Quantity):
        nframes (int, optional):
        nframes_limit (int, optional):
        dithers (list, optional):
        outdir (str, optional):
        outfile (str, optional):
        time_stamp (str, optional):
        debug (bool, optional):
    """

    # Compute number of files
    nfiles = nframes // nframes_limit
    nframes_left = nframes % nframes_limit
    logger.info("Creating {} files with {} synthetic exposures and adding {} frames to an additional file.".format(nfiles, nframes_limit, nframes_left))


    # Add a time stamp to file name, if not None
    outdir, outfile = os.path.split(outfile)
    now = datetime.now()
    time_str = now.strftime('%Y%m%d_%H%M%S')
    if time_stamp == 'end':
        outfile = outfile.replace('.fits', '_{}.fits'.format(time_str))
    elif time_stamp == 'start':
        outfile =  '{}_{}'.format(time_str, outfile)
    elif time_stamp is None:
        pass
    outfile = os.path.join(outdir, outfile)


    # Initialize fits header
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((nframes_limit,) + detector.shape)
    hdu.header.set('DIT', DIT.value, DIT.unit)
    hdu.header.set('DATE', str(datetime.now()))
    if 'cards' in  kwargs:
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
            card = "HIERARCH SPECKLEPY {} {}".format(object_name.upper(), key.upper())
            if isinstance(dict[key], u.Quantity):
                hdu.header.set(card, "{:.3e}".format(dict[key].value), dict[key].unit) # Appending the unit of a u.Quantity to the comment
            elif isinstance(dict[key], str):
                hdu.header.set(card, os.path.basename(dict[key])) # Suppress (long) relative paths
            elif isinstance(dict[key], tuple):
                _tuple = tuple([x.value for x in dict[key]]) # Separating tuple unit from values
                hdu.header.set(card, str(_tuple), dict[key][0].unit)
            else:
                hdu.header.set(card, dict[key])


    # Write header to one or more files, depending on 'nframes' and 'nframes_limit'
    outfiles = []
    for n in range(nfiles):
        filename = outfile.replace('.fits', '_{}.fits'.format(n + 1))
        outfiles.append(filename)
        hdu.writeto(filename, overwrite=True)
    # Create file for the left over frames
    if nframes_left > 0:
        filename = outfile.replace('.fits', '_{}.fits'.format(nfiles + 1))
        outfiles.append(filename)
        hdu.data = np.zeros((nframes_left,) + detector.shape)
        hdu.writeto(filename, overwrite=True)


    # Initialize parameters for frame computation
    if ('readout_time' in kwargs):
        skip_frames = int( kwargs['readout_time'] / telescope.psf_timestep )
    else:
        skip_frames = 0


    # Computation of frames
    frame_counter = 0
    for outfile_index, outfile in enumerate(outfiles):
        with fits.open(outfile, mode='update') as hdulist:
            # Get a new field of view for each file to enable dithering between files
            if dithers is not None:
                dither = dithers[outfile_index]
            else:
                dither = None
            photon_rate_density = target.get_photon_rate_density(FoV=detector.FoV, resolution=telescope.psf_resolution, dither=dither)

            for index in range(hdulist[0].header['NAXIS3']):
                photon_rate = telescope.get_photon_rate(photon_rate_density, integration_time=DIT, debug=debug)
                counts = detector.get_counts(photon_rate=photon_rate, integration_time=DIT, photon_rate_resolution=target.resolution, debug=debug)
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
        parameterfile (str): File from which the objects are instantiated.

    Returns:
        objects (dict): Dict containing the parameter file section as a key word
            and the objects (or kwargs dict) as values.
    """

    # Check whether files exist
    if not os.path.isfile(parameterfile):
        raise FileNotFoundError("Parameter file {} not found!".format(parameterfile))

    # Prepare objects list
    objects = {}

    # Read parameter_file
    parser = ConfigParser(inline_comment_prefixes="#")
    parser.optionxform = str  # make option names case sensitive
    logger.info("Reading parameter file {}".format(parameterfile))
    parser.read(parameterfile)
    for section in parser.sections():
        kwargs = {}
        for key in parser[section]:
            value = parser[section][key]
            try:
                kwargs[key] = eval(value)
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
        elif section.lower() == 'kwargs':
            objects['kwargs'] = kwargs

    return objects
