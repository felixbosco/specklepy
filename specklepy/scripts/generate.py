import os
from tqdm import trange

from astropy.io.fits import Header
from astropy.units import Quantity

from specklepy.io import FileStream, Config
from specklepy.logging import logger
from specklepy.mock import Detector, Target, Telescope
from specklepy.utils.time import default_time_stamp


def generate_exposure(target, telescope, detector, exposure_time, number_files=1, number_frames=1, dithers=None,
                      out_file='exposure.fits', time_stamp=None, debug=False, **kwargs):
    """Generate mock exposures from target, telescope and detector objects.

    The function generate_exposure() is the central function of the mock module. It takes one instance each  of
    the classes Target, Telescope, and Detector and the discrete exposure time. Then, it creates `number_files` files
    with `number_frames` frame each.

    Args:
        target (specklepy.mock.Target):
            Target instance that will be 'observed'.
        telescope (specklepy.mock.Telescope):
            Telescope instance that will 'observe'.
        detector (specklepy.mock.Detector):
            Detector instance that will be 'exposed'.
        exposure_time (astropy.units.Quantity):
            Discrete integration time for each exposure.
        number_files (int, optional):
            Number of files that will be generated. Default is 1.
        number_frames (int, optional):
            Number of frames per created file. Default is 1.
        dithers (list, optional):
            List of dither positions, offsets from the Target instance center. Default is None.
        out_file (str, optional):
            Base name of the output file. If multiple files are requested, the frame_index will be appended. Default is
            'exposure.fits'.
        time_stamp (str, optional):
            Time stamp that can be added to outfile. Can be either 'start', 'end' or None to suppress adding the time
            stamp. Default is None.
        debug (bool, optional):
            Show debugging information. Default is False.
    """

    #
    logger.info(f"Generating {number_files} files with {number_frames} mock exposure frames")

    # Add a time stamp to file name, if not None
    out_dir, out_file = os.path.split(out_file)
    time_str = default_time_stamp()
    if time_stamp == 'end':
        out_file = out_file.replace('.fits', f"_{time_str}.fits")
    elif time_stamp == 'start':
        out_file = f"{time_str}_{out_file}"
    elif time_stamp is None:
        pass
    out_file = os.path.join(out_dir, out_file)

    # Initialize fits header
    logger.info("Building FITS header...")
    header = Header()
    header.set('FEXPTIME', exposure_time.value, exposure_time.unit)
    header.set('EXPTIME', exposure_time.value * number_frames, exposure_time.unit)
    # header.set('DATE', str(datetime.now()))
    header.set('DATE', default_time_stamp())
    if 'cards' in kwargs:
        for key in kwargs['cards']:
            header.set(key, kwargs['cards'][key])
    # Add object attributes to header information
    skip_attributes = {'target': ['data', 'star_table', 'photometric_system'],
                       'telescope': ['psf', 'psf_frame'],
                       'detector': ['shape', 'array']}
    for _object in [target, telescope, detector]:
        object_dict = _object.__dict__
        object_name = _object.__name__

        for key in object_dict:
            if key in skip_attributes[object_name] or object_dict[key] is None:
                continue
            card = f"HIERARCH SPECKLEPY {object_name.upper()} {key.upper()}"
            if isinstance(object_dict[key], Quantity):
                # Appending the unit of a Quantity to the comment
                header.set(card, f"{object_dict[key].value:.3e}", object_dict[key].unit)
            elif isinstance(object_dict[key], str):
                header.set(card, os.path.basename(object_dict[key]))  # Suppress (long) relative paths
            elif isinstance(object_dict[key], tuple):
                _tuple = tuple([x.value for x in object_dict[key]])  # Separating tuple unit from values
                header.set(card, str(_tuple), object_dict[key][0].unit)
            else:
                header.set(card, object_dict[key])

    # Initialize parameters for frame computation
    if 'readout_time' in kwargs:
        skip_frames = int(kwargs['readout_time'] / telescope.psf_timestep)
    else:
        skip_frames = 0

    # Iterate through requested number of files
    for out_file_index in range(number_files):
        if number_files > 1:
            out_file_name = out_file.replace('.fits', f"_{out_file_index + 1}.fits")
        else:
            out_file_name = out_file
        file_stream = FileStream(out_file_name)
        file_stream.initialize(shape=(number_frames,) + detector.shape, header=header)

        # Get a new field of view for each file to enable dithering between files
        if dithers is not None:
            try:
                dither = dithers[out_file_index]
            except IndexError:
                raise RuntimeError(f"Expected {number_files} dither positions but received only {len(dithers)}!")
        else:
            dither = None
        logger.info("Computing photon rate density...")
        photon_rate_density = target.get_photon_rate_density(field_of_view=detector.field_of_view,
                                                             resolution=telescope.psf_resolution,
                                                             dither=dither)

        # Generate individual frames
        logger.info(f"Create {number_frames} exposures")
        for frame_index in trange(number_frames, desc=f"Generating exposures for file {out_file_name}"):
            photon_rate = telescope.get_photon_rate(photon_rate_density, integration_time=exposure_time,
                                                    debug=debug)
            counts = detector.get_counts(photon_rate=photon_rate,
                                         integration_time=exposure_time,
                                         photon_rate_resolution=target.resolution, debug=debug)
            counts = counts.decompose()
            file_stream.update_frame(frame_index=frame_index, data=counts.value)

            try:
                telescope.psf_frame += skip_frames
            except TypeError:
                pass


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
    params = Config.read(parameter_file).params

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
