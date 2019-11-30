import os
import numpy as np
from astropy.io import fits
import astropy.units as u
from datetime import datetime

from specklepy.logging import logging


def generate_exposure(target, telescope, detector, DIT, nframes=1, outdir=None, outfile=None, time_stamp='end', nframes_limit=100, verbose=0, **kwargs):
    """Generate synthetic exposures from target, telescope and detector objects.

    The function generate_exposure() is the central function of the synthetic
    module. It takes one instance each  of the classes Target, Telescope, and
    Detector and the discrete exposure time DIT. Then, it creates a number of
    files depending on the number of requested frames ('nframes') and frame
    limit per file ('nframes_limit').
    To distribute the synthetic exposures to multiple files, for instance if the
    size of the individual file would become too large, just set
    'nframes_limit' to a smaller value (default is 100). The
    last file may contain empty frames.

    Args:
        target ():
        telescope ():
        detector ():
        DIT ():
        ...
    """

    # Compute number of files
    nfiles = nframes // nframes_limit
    nframes_left = nframes % nframes_limit
    logging.info("Creating {} files with {} synthetic exposures and adding {} frames to an additional file.".format(nfiles, nframes_limit, nframes_left))

    # Adapt file name
    if outfile is None:
        outfile = 'exposure.fits'
    if time_stamp == 'end':
        try:
            generic, ext = outfile.split('.')
            outfile = generic + '_' + _make_time_stamp() + '.' + ext
        except ValueError as e:
            path = outfile
            outfile = outfile.split('/')[-1]
            generic, ext = outfile.split('.')
            outfile = path.replace(outfile, generic + '_' + _make_time_stamp() + '.' + ext)
    elif time_stamp == 'start':
        outfile =  _make_time_stamp() + '_' + outfile
    elif time_stamp is None:
        pass
    outfile = os.path.join(outdir, outfile)


    # Initialize fits header
    hdu = fits.PrimaryHDU()
    # hdu.header.set('NAXIS', 2)
    # hdu.header.set('NAXIS1', detector.shape[0])
    # hdu.header.set('NAXIS2', detector.shape[1])
    # if nframes > 1:
    #     # In case of multiple frames, update 'NAXIS'
    #     hdu.header.set('NAXIS', 3, 'number of array dimensions')
    #     hdu.header.set('NAXIS3', nframes)
    #     hdu.data = np.zeros( (nframes, detector.shape[0], detector.shape[1]) )
    # else:
    #     hdu.data = np.zeros(detector.shape)
    hdu.data = np.zeros((nframes_limit,) + detector.shape)
    hdu.header.set('DIT', DIT.value, DIT.unit)
    _add_attributes_to_header(hdu, target, skip_attributes=['shape', 'data', 'stars'], object_name='TARGET')
    _add_attributes_to_header(hdu, telescope, skip_attributes=['psf'], object_name='TELESCOP')
    _add_attributes_to_header(hdu, detector, skip_attributes=['shape', 'array'], object_name='DETECTOR')
    hdu.header.set('DATE', str(datetime.now()))


    # Write header to one or more files, depending on 'nframes' and 'nframes_limit'
    # if nframes <= nframes_limit:
    #     multiple_files = False
    #     print("Writing file {}.".format(outfile))
    #     hdu.writeto(outfile, overwrite=True)
    # else:
    #     multiple_files = True
    #     number_full_files = nframes // nframes_limit
    #     number_leftover_frames = nframes % nframes_limit
    #     if number_leftover_frames != 0:
    #         print("Writing {} files, where the last file contains only {} valid frames.".format(number_full_files + 1, number_leftover_frames))
    #
    #         # Writing files with the maximum number of frames
    #         for i in range(number_full_files):
    #             hdu.header.set('NAXIS3', nframes_limit)
    #             hdu.writeto(_make_outfile(outfile, i, add_index=multiple_files), overwrite=True)
    #
    #         # The last file shall contain only fewer frames
    #         hdu.header.set('NAXIS3', number_leftover_frames)
    #         hdu.writeto(_make_outfile(outfile, i+1, add_index=multiple_files), overwrite=True)
    #     else:
    #         print("Writing {} files.".format(number_full_files))
    #
    #         # Writing files with the maximum number of frames
    #         for i in range(number_full_files):
    #             hdu.header.set('NAXIS3', nframes_limit)
    #             hdu.writeto(_make_outfile(outfile, i, add_index=multiple_files), overwrite=True)
    outfiles = []
    for n in range(nfiles):
        filename = outfile.replace('.fits', '_{}.fits'.format(n + 1))
        outfiles.append(filename)
        hdu.writeto(filename, overwrite=True)
    if nframes_left > 0:
        # Create file for the left over frames
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
    for outfile in outfiles:
        with fits.open(outfile, mode='update') as hdulist:
            for index in range(hdulist[0].header['NAXIS3']):
                imaged = telescope(target.data, target.resolution, integration_time=DIT, verbose=verbose)
                detected = detector(photon_rate_density_array=imaged, integration_time=DIT, target_FoV=target.FoV)
                detected = detected.decompose()
                hdulist[0].data[index] = detected.value

                try:
                    telescope.psf_plane += skip_frames
                except TypeError:
                    pass

                frame_counter += 1
                print("\rExposure {:4}/{:4}".format(frame_counter, nframes), end='')
    print("")

    # for dt in range(nframes):
    #     print("\rExposure {:4}/{:4}".format(dt+1, nframes), end='')
    #     imaged = telescope(target.data, target.resolution, integration_time=DIT, verbose=verbose)
    #     detected = detector(photon_rate_density_array=imaged, integration_time=DIT, target_FoV=target.FoV)
    #     detected = detected.decompose()
    #     # Write file
    #     with fits.open(_make_outfile(outfile, dt // nframes_limit, add_index=multiple_files), mode='update') as hdulist:
    #         if nframes == 1:
    #             hdulist[0].data = detected.value
    #         else:
    #             if multiple_files:
    #                 hdulist[0].data[dt % nframes_limit] = detected.value
    #             else:
    #                 hdulist[0].data[dt] = detected.value
    #         hdulist.flush()
    #     # Skip psf frames, to account for time between two readouts
    #     try:
    #         telescope.psf_plane += skip_frames
    #     except TypeError:
    #         pass
    # print("")



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

def _add_attributes_to_header(hdu_object, object, skip_attributes=[], prefix='HIERARCH SPECKLEPY ', object_name='New object'):
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

def _make_outfile(outfile, index, add_index):
    if add_index:
        generic, extension = outfile.split('.')
        return "{}_{}.{}".format(generic, index, extension)
    else:
        return outfile
