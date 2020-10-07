import os
from datetime import datetime
import numpy as np
from tqdm import trange

from astropy.io import fits
from astropy.units import Quantity

from specklepy.io import config
from specklepy.logging import logger
from specklepy.synthetic.target import Target
from specklepy.synthetic.telescope import Telescope
from specklepy.synthetic.detector import Detector


def generate_exposure(target, telescope, detector, exposure_time, n_frames=1, n_frames_limit=100, dithers=None,
                      outfile='exposure.fits', time_stamp=None, debug=False, **kwargs):
    """Generate synthetic exposures from target, telescope and detector objects.

    The function generate_exposure() is the central function of the synthetic module. It takes one instance each  of
    the classes Target, Telescope, and Detector and the discrete exposure time. Then, it creates a number of files
    depending on the number of requested frames ('n_frames') and frame limit per file ('n_frames_limit'). To distribute
    the synthetic exposures to multiple files, for instance if the size of the individual file would become too large,
    just set 'n_frames_limit' to a smaller value (default is 100).

    Args:
        target (specklepy.synthetic.Target):
            Target instance that will be 'observed'.
        telescope (specklepy.synthetic.Telescope):
            Telescope instance that will 'observe'.
        detector (specklepy.synthetic.Detector):
            Detector instance that will be 'exposed'.
        exposure_time (astropy.units.Quantity):
            Discrete integration time for each exposure.
        n_frames (int, optional):
            Number of frames that will be generated. Default is 1.
        n_frames_limit (int, optional):
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
    n_files = n_frames // n_frames_limit
    n_frames_left = n_frames % n_frames_limit
    if n_frames_left:
        logger.info(f"Creating {n_files} files with {n_frames_limit} synthetic exposures and adding {n_frames_left} "
                    f"frames to an additional file")
    else:
        logger.info(f"Creating {n_files} files with {n_frames_limit} synthetic exposures")

    # Add a time stamp to file name, if not None
    out_dir, outfile = os.path.split(outfile)
    now = datetime.now()
    time_str = now.strftime('%Y%m%d_%H%M%S')
    if time_stamp == 'end':
        outfile = outfile.replace('.fits', f"_{time_str}.fits")
    elif time_stamp == 'start':
        outfile = f"{time_str}_{outfile}"
    elif time_stamp is None:
        pass
    outfile = os.path.join(out_dir, outfile)

    # Initialize fits header
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((n_frames_limit,) + detector.shape)
    hdu.header.set('EXPTIME', exposure_time.value, exposure_time.unit)
    hdu.header.set('DATE', str(datetime.now()))
    if 'cards' in kwargs:
        for key in kwargs['cards']:
            hdu.header.set(key, kwargs['cards'][key])
    # Add object attributes to header information
    skip_attributes = {'target': ['data', 'stars'],
                       'telescope': ['psf', 'psf_frame'],
                       'detector': ['shape', 'array']}
    for _object in [target, telescope, detector]:
        object_dict = _object.__dict__
        object_name = _object.__name__

        for key in object_dict:
            if key in skip_attributes[object_name]:
                continue
            card = f"HIERARCH SPECKLEPY {object_name.upper()} {key.upper()}"
            if isinstance(object_dict[key], Quantity):
                # Appending the unit of a Quantity to the comment
                hdu.header.set(card, f"{object_dict[key].value:.3e}", object_dict[key].unit)
            elif isinstance(object_dict[key], str):
                hdu.header.set(card, os.path.basename(object_dict[key]))  # Suppress (long) relative paths
            elif isinstance(object_dict[key], tuple):
                _tuple = tuple([x.value for x in object_dict[key]])  # Separating tuple unit from values
                hdu.header.set(card, str(_tuple), object_dict[key][0].unit)
            else:
                hdu.header.set(card, object_dict[key])

    # Write header to one or more files, depending on 'n_frames' and 'n_frames_limit'
    outfiles = []
    for n in range(n_files):
        filename = outfile.replace('.fits', f"_{n + 1}.fits")
        outfiles.append(filename)
        hdu.writeto(filename, overwrite=True)
    # Create file for the left over frames
    if n_frames_left > 0:
        filename = outfile.replace('.fits', f"_{n_files + 1}.fits")
        outfiles.append(filename)
        hdu.data = np.zeros((n_frames_left,) + detector.shape)
        hdu.writeto(filename, overwrite=True)

    # Initialize parameters for frame computation
    if 'readout_time' in kwargs:
        skip_frames = int(kwargs['readout_time'] / telescope.psf_timestep)
    else:
        skip_frames = 0

    # Computation of frames
    frame_counter = 0
    for outfile_index, outfile in enumerate(outfiles):
        with fits.open(outfile, mode='update') as hdu_list:

            # Get a new field of view for each file to enable dithering between files
            if dithers is not None:
                try:
                    dither = dithers[outfile_index]
                except IndexError:
                    raise RuntimeError(f"Expected {len(outfiles)} dither positions but received only {len(dithers)}!")
            else:
                dither = None
            photon_rate_density = target.get_photon_rate_density(field_of_view=detector.field_of_view,
                                                                 resolution=telescope.psf_resolution,
                                                                 dither=dither)

            # Generate individual frames
            for index in trange(hdu_list[0].header['NAXIS3'], desc=f"Generating exposures for file {outfile}"):
                photon_rate = telescope.get_photon_rate(photon_rate_density, integration_time=exposure_time,
                                                        debug=debug)
                counts = detector.get_counts(photon_rate=photon_rate,
                                             integration_time=exposure_time,
                                             photon_rate_resolution=target.resolution, debug=debug)
                counts = counts.decompose()
                hdu_list[0].data[index] = counts.value

                try:
                    telescope.psf_frame += skip_frames
                except TypeError:
                    pass

                frame_counter += 1

            # Update header entry DATE
            hdu_list[0].header.set('DATE', str(datetime.now()))


def get_objects(parameter_file, debug=False):
    """Get objects from parameter file.

    Args:
        parameter_file (str):
            File from which the objects are instantiated.
        debug (bool, optional):
            Show debugging information.

    Returns:
        target (Target object):
        telescope (Telescope object):
        detector (Detector object):
        parameters (dict):
            Dictionary containing the parameters parsed toward generate exposure beyond the three objects.
    """

    # Set logger level
    if debug:
        logger.setLevel('DEBUG')

    # Check whether files exist
    if not os.path.isfile(parameter_file):
        raise FileNotFoundError(f"Parameter file {parameter_file} not found!")

    # Read parameter file
    params = config.read(parameter_file)

    # Create objects from the parameters
    target = Target(**params['TARGET'])
    logger.debug(f"Initialized Target instance:\n{target}")
    telescope = Telescope(**params['TELESCOPE'])
    logger.debug(f"Initialized Telescope instance:\n{telescope}")
    detector = Detector(**params['DETECTOR'])
    logger.debug(f"Initialized Detector instance:\n{detector}")

    # Extract and interpret other kez word parameters
    if 'KWARGS' in params:
        parameters = params['KWARGS']
    elif 'PARAMETERS' in params:
        parameters = params['PARAMETERS']
    else:
        parameters = {}

    # Interpret str-type key word arguments
    for key in parameters.keys():
        if isinstance(parameters[key], str):
            try:
                parameters[key] = eval(parameters[key])
                logger.debug(f"Kwarg {key} evaluated as {parameters[key]} ({type(parameters[key])})")
            except SyntaxError:
                parameters[key] = Quantity(parameters[key])
                logger.debug(f"Kwarg {key} evaluated as {parameters[key]} ({type(parameters[key])})")
            except NameError:
                logger.debug(f"Kwarg {key} not evaluated {parameters[key]} ({type(parameters[key])})")

    return target, telescope, detector, parameters
