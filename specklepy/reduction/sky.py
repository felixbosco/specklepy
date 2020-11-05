from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from tqdm import trange
import os

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats

from photutils import make_source_mask

from specklepy.io.outfile import Outfile
from specklepy.logging import logger
from specklepy.exceptions import SpecklepyTypeError
from specklepy.plotting.plot import save_figure


def subtract_sky_background(in_files, out_files=None, method='scalar', source='sky', mask_sources=False, file_path=None,
                            tmp_dir=None, show=False, debug=False):

    """Estimate and subtract the sky background via different methods and sources.

    TODO: Implement sky subtraction from image

    Args:
        in_files (specklepy.FileArchive):
            File archive storing the information of all the files in the reduction.
        out_files (list):
            List of files to apply the sky subtraction to. If left empty, the list stored in the `in_files` FileArchive
            is used.
        method (str, optional):
            Switch between a scalar (`scalar`) background value or a 2D image (`images`).
        source (str, optional):
            Source for estimating the background from. Can be either `sky` to measure from dedicated sky frames or
            `science` to use the science frames themselves. Typically, these frames have a high number of sources, so
            `mask_sources` should be switched on.
        mask_sources (bool, optional):
            In empty reference fields, this masking option should stay at `False`, since source masking is not well
            tested. However, masking sources yields a more precise result.
        file_path (str, optional):
            Path to the files, listed in `in_files`.
        tmp_dir (str, optional):
            Directory to which temporary results and QA data is stored.
        show (bool, optional):
            Show plots of sky estimates for each sequence. They will be created and stored regardless of this choice.
        debug (bool, optional):
            Show debugging information.
    """

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')

    # Apply fall back values
    if method is None:
        method = 'scalar'
    logger.info(f"Sky background subtraction method: {method}")
    if source is None:
        source = 'sky'
    logger.info(f"Sky background subtraction source: {source}")
    if out_files is None:
        out_files = in_files.product_files
    if out_files is None:
        logger.warning(f"Output files are not declared in subtract_sky_background!")

    # Identify the observing sequences
    sequences = in_files.identify_sequences(source=source)

    # Start the background estimates
    if method == 'scalar':

        # Iterate through observing sequences
        for s, sequence in enumerate(sequences):

            logger.info(f"Starting observing sequence {s} :: Object {sequence.object} :: Setup {sequence.setup}")

            # Compute weights based on the time offset to the individual sky observations
            weights = sequence.compute_weights()

            # Start extracting sky fluxes
            sky_bkg = np.zeros(sequence.n_sky)
            sky_bkg_std = np.zeros(sequence.n_sky)
            for i in trange(sequence.n_sky, desc='Estimate sky background from cube'):
                file = sequence.sky_files[i]
                bkg, d_bkg = estimate_sky_background(file, method=method, mask_sources=mask_sources, path=file_path)
                sky_bkg[i] = bkg
                sky_bkg_std[i] = d_bkg
            logger.debug(f"Shapes:\nF: {sky_bkg.shape}\ndF: {sky_bkg_std.shape}")

            # Compute weighted sky background for each science file
            weighted_sky_bkg = np.dot(weights, sky_bkg)
            weighted_sky_bkg_var = np.dot(np.square(weights), np.square(sky_bkg_std))

            # Store sky background estimates
            sky_bkg_table = Table(data=[sequence.sky_files, weighted_sky_bkg, weighted_sky_bkg_var],
                                  names=['FILE', 'BKG', 'VAR'])
            sky_bkg_table_name = f"sky_bkg_{sequence.object}_{sequence.setup}.fits"
            sky_bkg_table.write(os.path.join(tmp_dir, sky_bkg_table_name), overwrite=True)

            # Plot sky flux estimates
            for i, file in enumerate(sequence.sky_files):
                plt.text(sequence.sky_time_stamps[i], sky_bkg[i], file, rotation=90, alpha=.5)
            for i, file in enumerate(sequence.science_files):
                plt.text(sequence.science_time_stamps[i], weighted_sky_bkg[i], file, rotation=90,
                         alpha=.66)
            plt.errorbar(x=sequence.sky_time_stamps, y=sky_bkg, yerr=sky_bkg_std,
                         fmt='None', ecolor='tab:blue', alpha=.5)
            plt.plot(sequence.sky_time_stamps, sky_bkg, 'D', label='Sky', c='tab:blue')
            plt.errorbar(x=sequence.science_time_stamps, y=weighted_sky_bkg, yerr=np.sqrt(weighted_sky_bkg_var),
                         fmt='None', ecolor='tab:orange', alpha=.66)
            plt.plot(sequence.science_time_stamps, weighted_sky_bkg, 'D', label='Science', c='tab:orange')
            plt.xlabel('Time (s)')
            plt.ylabel('Flux (counts)')
            plt.legend()
            save_figure(os.path.join(tmp_dir, sky_bkg_table_name.replace('.fits', '.png')))
            if show:
                plt.show()
            plt.close()

            # Subtract sky from product files
            for i, science_file in enumerate(sequence.science_files):
                for out_file in out_files:
                    if science_file in out_file:
                        science_file = out_file
                logger.info(f"Applying sky background subtraction on file {science_file}")
                with fits.open(science_file, mode='update') as hdu_list:
                    hdu_list[0].data = hdu_list[0].data.astype(float) - weighted_sky_bkg[i]
                    if 'VAR' in hdu_list:
                        hdu_list['VAR'].data = hdu_list['VAR'].data + weighted_sky_bkg_var[i]
                    else:
                        # Construct new HDU
                        shape = np.array(hdu_list[0].data.shape)[[-2, -1]]
                        data = np.full(shape=shape, fill_value=weighted_sky_bkg_var[i])
                        hdu = fits.ImageHDU(data=data, name='VAR')
                        hdu_list.append(hdu)
                    hdu_list[0].header.set('SKYCORR', str(datetime.now()))
                    hdu_list[0].header.set('SKYBKG', weighted_sky_bkg[i], "Sky background")
                    hdu_list[0].header.set('SKYVAR', weighted_sky_bkg_var[i], "Sky background variance")
                    hdu_list.flush()

    elif method in ['image', 'frame']:
        raise NotImplementedError("Sky subtraction in image mode is not implemented yet!")

    else:
        raise ValueError(f"Sky subtraction method {method} is not understood!")


def subtract_scalar_background(files, params, prefix=None, debug=False):
    """Estimate and subtract a scalar background."""

    if not isinstance(files, (list, np.ndarray)):
        raise SpecklepyTypeError('substract_scalar_background', argtype=type(files), argname='files', expected='list')
    else:
        if len(files) == 0:
            raise RuntimeError("Sky subtraction received an empty list of files!")

    logger.info("Estimating scalar background and subtract...")
    for file_index, file in enumerate(files):

        image, header = fits.getdata(os.path.join(params.paths.filePath, file), header=True)

        # Update header for and initialize the outfile
        header.set('PIPELINE', 'SPECKLEPY')
        header.set('SKYCORR', str(datetime.now()))
        corrected_file = prefix + file
        corrected_file = os.path.join(params.paths.filePath, corrected_file)
        outfile = Outfile(filename=corrected_file, header=header, shape=image.shape)

        # Estimate scalar background and uncertainties, and subtract
        if image.ndim == 2:
            mean, median, std = sigma_clipped_stats(image, sigma=params.sky.backgroundSigmaClip)
            outfile.data = image - mean
            image_var = np.ones(image.shape) * np.square(std)
            outfile.new_extension(name='VAR', data=image_var)
        elif image.ndim == 3:
            means, medians, stds = sigma_clipped_stats(image, sigma=params.sky.backgroundSigmaClip, axis=(1, 2))
            logger.info(f"Sigma clipped stats:\t{np.mean(means):.2f} +- {np.mean(stds):.2f}")
            outfile.new_extension(name='VAR', data=np.zeros(image.shape))
            tmp_frame = np.ones(image[0].shape)

            for frame_index, frame in enumerate(image):
                print(f"\r\tUpdating frame {frame_index+1:3}...", end='')
                outfile.update_frame(frame_index=frame_index, data=np.subtract(frame, means[frame_index]))
                outfile.update_frame(frame_index=frame_index, data=tmp_frame * np.square(stds[frame_index]),
                                     extension='VAR')
            print()
        else:
            raise RuntimeError(f"Images are supposed to have 2 or 3 dimensions but this one has {image.ndim}!")

    logger.info("Scalar background subtraction complete!")


def estimate_sky_background(data, method='scalar', mask_sources=True, path=None):

    """Estimate a scalar or image sky background with uncertainties.

    Args:
        data (np.array or str):
            Image or data cube with sky observations used to extract the sky background mean and uncertainty. Str type
            input is interpreted as a file name to read the data from.
        method (str, optional):
            Can be `scalar` or `image` for setting the shape of the output.
        mask_sources (bool, optional):
            Creates a source mask to exclude sources from the measurement. This should not be set if working in `image`
            mode.
        path (str, optional):
            Path to the data files. This is used only if data is provided as a file name.

    Returns:
        mean (float or np.array):
            Mean sky background as a scalar or image, depending on the method parameter.
        std (float or np.array):
            Uncertainty on the sky background estimate as a scalar or image, depending on the method parameter.
    """

    # Handle str type data
    if isinstance(data, str):
        file = data
        if path is not None:
            file = os.path.join(path, file)
        data = fits.getdata(filename=file)

    # Create source mask
    if mask_sources:
        if data.ndim == 3:
            mask = make_source_mask(np.sum(data, axis=0), nsigma=2, npixels=5, dilate_size=11)
            mask = np.repeat(np.expand_dims(mask, 0), data.shape[0], axis=0)
        else:
            mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
    else:
        mask = None

    # Derive statistics
    if method == 'scalar':
        mean, _, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    else:
        raise NotImplementedError(f"Method {method} is not implemented yet for sky background estimation!")

    return mean, std
