import os
import numpy as np
from astropy.io import fits
import astropy.units as u
from datetime import datetime

from specklepy.logging import logging


def generate_exposure(target, telescope, detector, DIT, nframes=1, nframes_limit=100, outdir=None, outfile='exposure.fits', time_stamp='end', verbose=False, **kwargs):
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
        outdir (str, optional):
        outfile (str, optional):
        time_stamp (str, optional):
        verbose (bool, optional):
    """

    # Compute number of files
    nfiles = nframes // nframes_limit
    nframes_left = nframes % nframes_limit
    logging.info("Creating {} files with {} synthetic exposures and adding {} frames to an additional file.".format(nfiles, nframes_limit, nframes_left))


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
    # Add object attributes to header information
    skip_attributes = {'target': ['shape', 'data', 'stars'],
                       'telescope': ['psf'],
                       'detector': ['shape', 'array']}
    for object in [target, telescope, detector]:
        dict = object.__dict__
        object_name = object.__name__

        for key in dict:
            if key in skip_attributes[object_name]:
                continue
            card = "HIERARCH SPECKLEPY {} {}".format(object_name.upper(), key.upper())
            if isinstance(dict[key], u.Quantity):
                hdu.header.set(card, dict[key].value, dict[key].unit) # Appending the unit of a u.Quantity to the comment
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
                print("\rSaving exposure {:4}/{:4} to file {}...".format(frame_counter, nframes, outfile), end='')
    print("")
