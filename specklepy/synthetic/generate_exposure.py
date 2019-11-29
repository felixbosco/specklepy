# Function generate_exposure()
# The function makes use of several helper functions, which can be found below the function definition.



# Dependencies
import os
import numpy as np
from astropy.io import fits
import astropy.units as u
from datetime import datetime



# Main definition
def generate_exposure(target, telescope, detector, DIT, number_frames=1, outdir=None, filename=None, time_stamp='end', maximum_number_frames_per_file=100, verbose=0, **kwargs):
    """
    The function generate_exposure() is the central function of the VEGAPy
    package.
    It takes three quantities carrrying information on the vegapy.Target,
    vegapy.Telescope, and vegapy.Detector. It generates 'number_frames'
    exposures of integration time 'DIT' and writes them to a fits file
    'filename'. The filename automatically obtains a time stamp, as long as the
    argument is set to 'start' or to the default 'end'.
    To distribute the virtual exposures to multiple files, for instance if the
    size of the file would become too large, just set
    'maximum_number_frames_per_file' to a smaller value (default is 100). The
    last file may contain empty frames.
    """

    # Adapt file name
    if filename is None:
        filename = 'exposure.fits'
    if time_stamp == 'end':
        try:
            generic, ext = filename.split('.')
            filename = generic + '_' + _make_time_stamp() + '.' + ext
        except ValueError as e:
            path = filename
            filename = filename.split('/')[-1]
            generic, ext = filename.split('.')
            filename = path.replace(filename, generic + '_' + _make_time_stamp() + '.' + ext)
    elif time_stamp == 'start':
        filename =  _make_time_stamp() + '_' + filename
    elif time_stamp is None:
        pass
    filename = os.path.join(outdir, filename)


    # Initialize fits header
    hdu = fits.PrimaryHDU()
    hdu.header.set('NAXIS', 2)
    hdu.header.set('NAXIS1', detector.shape[0])
    hdu.header.set('NAXIS2', detector.shape[1])
    if number_frames > 1:
        # In case of multiple frames, update 'NAXIS'
        hdu.header.set('NAXIS', 3, 'number of array dimensions')
        hdu.header.set('NAXIS3', number_frames)
        hdu.data = np.zeros( (number_frames, detector.shape[0], detector.shape[1]) )
    else:
        hdu.data = np.zeros(detector.shape)
    hdu.header.set('DIT', DIT.value, DIT.unit)
    _add_attributes_to_header(hdu, target, skip_attributes=['shape', 'data', 'stars'], object_name='TARGET')
    _add_attributes_to_header(hdu, telescope, skip_attributes=['psf'], object_name='TELESCOP')
    _add_attributes_to_header(hdu, detector, skip_attributes=['shape', 'array'], object_name='DETECTOR')
    hdu.header.set('DATE', str(datetime.now()))


    # Write header to one or more files, depending on 'number_frames' and 'maximum_number_frames_per_file'
    if number_frames <= maximum_number_frames_per_file:
        multiple_files = False
        print("Writing file {}.".format(filename))
        hdu.writeto(filename, overwrite=True)
    else:
        multiple_files = True
        number_full_files = number_frames // maximum_number_frames_per_file
        number_leftover_frames = number_frames % maximum_number_frames_per_file
        if number_leftover_frames != 0:
            print("Writing {} files, where the last file contains only {} valid frames.".format(number_full_files + 1, number_leftover_frames))

            # Writing files with the maximum number of frames
            for i in range(number_full_files):
                hdu.header.set('NAXIS3', maximum_number_frames_per_file)
                hdu.writeto(_make_filename(filename, i, add_index=multiple_files), overwrite=True)

            # The last file shall contain only fewer frames
            hdu.header.set('NAXIS3', number_leftover_frames)
            hdu.writeto(_make_filename(filename, i+1, add_index=multiple_files), overwrite=True)
        else:
            print("Writing {} files.".format(number_full_files))

            # Writing files with the maximum number of frames
            for i in range(number_full_files):
                hdu.header.set('NAXIS3', maximum_number_frames_per_file)
                hdu.writeto(_make_filename(filename, i, add_index=multiple_files), overwrite=True)




    # Initialize parameters for frame computation
    if ('readout_time' in kwargs):
        skip_frames = int( kwargs['readout_time'] / telescope.psf_timestep )
    else:
        skip_frames = 0


    # Computation of frames
    for dt in range(number_frames):
        print("\rExposure {:4}/{:4}".format(dt+1, number_frames), end='')
        imaged = telescope(target.data, target.resolution, integration_time=DIT, verbose=verbose)
        detected = detector(photon_rate_density_array=imaged, integration_time=DIT, target_FoV=target.FoV)
        detected = detected.decompose()
        # Write file
        with fits.open(_make_filename(filename, dt // maximum_number_frames_per_file, add_index=multiple_files), mode='update') as hdulist:
            if number_frames == 1:
                hdulist[0].data = detected.value
            else:
                if multiple_files:
                    hdulist[0].data[dt % maximum_number_frames_per_file] = detected.value
                else:
                    hdulist[0].data[dt] = detected.value
            hdulist.flush()
        # Skip psf frames, to account for time between two readouts
        try:
            telescope.psf_plane += skip_frames
        except TypeError:
            pass
    print("")



# Helper functions
def _make_time_stamp():
    """
    The helper function _make_time_stamp() returns a string:
    'YYYYMMDD_HHMMSS'.
    """
    tmp = str(datetime.now())
    tmp = tmp.split('.')[0]
    tmp = tmp.replace(' ', '_')
    tmp = tmp.replace('-', '')
    tmp = tmp.replace(':', '')
    return tmp

def _add_attributes_to_header(hdu_object, object, skip_attributes=[], prefix='HIERARCH VEGAPY ', object_name='New object'):
    """
    The helper function _add_attributes_to_header() formats the attributes of
    the argument object into appropriate FITS header cards.
    For distinguishing the attributes of different objects, a headline card is
    created for the given object, followed by HIERARCH cards with the attributes
    as long as prefix is set to 'HIERARCH '.
    """
    dict = object.__dict__
    #hdu_object.header.set(object_name, '')
    for key in dict:
        # Ability to skip for instance arrays
        if key in skip_attributes:
            continue
        # Appending the unit of a u.Quantity to the comment
        if isinstance(dict[key], u.Quantity):
            hdu_object.header.set(prefix + object_name + ' ' + key, dict[key].value, dict[key].unit)
        # Suppress (long) relative paths
        elif isinstance(dict[key], str):
            if len(dict[key]) > 20:
                path, file = os.path.split(dict[key])
                hdu_object.header.set(prefix + object_name + ' ' + key, file)
            else:
                hdu_object.header.set(prefix + object_name + ' ' + key, dict[key])
        # Separating tuple attributes into two header cards
        elif isinstance(dict[key], tuple):
            hdu_object.header.set(prefix + object_name + ' ' + key + '[0]', dict[key][0].value, dict[key][0].unit)
            hdu_object.header.set(prefix + object_name + ' ' + key + '[1]', dict[key][1].value, dict[key][1].unit)
        # Add all other types
        else:
            hdu_object.header.set(prefix + object_name + ' ' + key, dict[key])

def _make_filename(filename, index, add_index):
    if add_index:
        generic, extension = filename.split('.')
        return "{}_{}.{}".format(generic, index, extension)
    else:
        return filename
