import numpy as np
import os
import sys

from astropy.io import fits

from specklepy.core import alignment
from specklepy.core.sourceextraction import extract_sources
from specklepy.core.specklecube import SpeckleCube
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.fits import get_data
from specklepy.io.reconstructionfile import ReconstructionFile
from specklepy.logging import logger
from specklepy.utils.array import frame_shape
from specklepy.utils.box import Box


class Reconstruction(object):

    """Base class for image reconstructions.

    This class stores the principle parameters of the reconstructed image. In the future, the functions `ssa` and
    `holography` may become child classes of this one.
    """

    supported_modes = ['full', 'same', 'valid']
    variance_extension = 'VAR'
    mask_extension = 'MASK'
    fill_value = 0

    def __init__(self, in_files, mode='same', reference_file=None, out_file=None, in_dir=None, tmp_dir=None,
                 integration_method='collapse', variance_extension=None, box_indexes=None, custom_mask=None,
                 debug=False):
        """Create a Reconstruction instance.

        Args:
            in_files (list):
                List of input data cubes.
            mode (str, optional):
                Reconstruction mode, defines the final image size and can be `full`, `same` and `valid`. The final image
                sizes is derived as follows:
                - `full`:
                    The reconstruction image covers every patch of the sky that is covered by at least one frame in the
                    input data.
                - `same`:
                    The reconstruction image covers the same field of view as the image in the reference file.
                - `valid`:
                    The reconstruction image covers only that field that is covered by all images in the input files.
            reference_file (int or str, optional):
                The index in the `in_files` list or the name of the image serving as reference in 'same' mode.
            out_file (str, optional):
                Name of an output file to store the reconstructed image in.
            in_dir (str, optional):
                Path to the `in_files`.
            tmp_dir (str, optional):
                Path to the directory for storing temporary products.
            integration_method (str, optional):
                Method for integrating speckle cubes. Can also be set in the `create_long_exposures` method.
            variance_extension (str, optional):
                Name of the variance extension of the FITS files. The default is `'VAR'`.
            box_indexes (list, optional):
                List of indexes defining a box for the SSA integration.
            custom_mask (str or array-like, optional):
                Path to a custom mask or the mask itself, which will be applied to all frames prior to integration.
            debug (bool, optional):
                Show debugging information.
        """

        # Check input parameter types
        if not isinstance(in_files, (list, np.ndarray)):
            raise SpecklepyTypeError('Reconstruction', 'in_files', type(in_files), 'list')
        if not isinstance(mode, str):
            raise SpecklepyTypeError('Reconstruction', 'mode', type(mode), 'str')
        if out_file is not None and not isinstance(out_file, str):
            raise SpecklepyTypeError('Reconstruction', 'out_file', type(out_file), 'str')

        # Check input parameter values
        if mode not in self.supported_modes:
            raise SpecklepyValueError('Reconstruction', 'mode', mode, f"in {self.supported_modes}")

        # Store input data
        self.in_files = in_files
        self.mode = mode
        self.out_file = out_file if out_file is not None else 'reconstruction.fits'
        self.reference_index = self.identify_reference_file(reference_file)
        self.in_dir = in_dir if in_dir is not None else ''
        self.tmp_dir = tmp_dir if tmp_dir is not None else ''
        self.integration_method = integration_method
        if variance_extension is not None:
            self.variance_extension = variance_extension  # if var_ext is not None else 'VAR'
        self.box = Box(box_indexes) if box_indexes is not None else None
        self.debug = debug

        # Set custom mask
        self.custom_mask_file = None
        self._custom_mask = None
        if isinstance(custom_mask, str):
            self.custom_mask_file = custom_mask
        elif isinstance(custom_mask, (bool, np.ndarray)):
            self._custom_mask = custom_mask

        # Derive shape of individual input frames
        example_data = get_data(self.in_files[0], path=in_dir)
        self.frame_shape = frame_shape(example_data)

        # Initialize secondary attributes
        self.image = None
        self.image_var = None
        self.long_exp_files = None
        self.shifts = None
        self.pad_vectors = None
        self.reference_pad_vector = None

        # Initialize output file and create an extension for the variance
        self.out_file = ReconstructionFile(files=self.in_files, filename=self.out_file,
                                           in_dir=in_dir, cards={"RECONSTRUCTION": "SSA"})
        if self.variance_extension is not None:
            self.out_file.new_extension(name=self.variance_extension)

    @property
    def single_cube_mode(self):
        return len(self.in_files) == 1

    @property
    def reference_file(self):
        return self.in_files[self.reference_index]

    @property
    def custom_mask(self):
        if self._custom_mask is not None:
            return self._custom_mask
        elif self.custom_mask_file is not None:
            return get_data(self.custom_mask_file, dtype=bool)
        else:
            return False

    def assert_dirs(self, error=True):
        for name, dir in zip(['in_dir', 'tmp_dir'], [self.in_dir, self.tmp_dir]):
            # Check for None-types
            if dir is None:
                if error:
                    raise ValueError(f"{name}-directory is not provided!")
                else:
                    continue

            # Create directory if not existing yet
            if not os.path.isdir(dir):
                logger.info(f"Creating {name} {dir!r}")
                os.makedirs(dir)

    def identify_reference_file(self, reference_file=None):

        # Initialize reference
        reference_index = None

        # Interpret reference image
        if reference_file is None:
            reference_index = 0
        elif isinstance(reference_file, int):
            reference_index = reference_file
        elif isinstance(reference_file, str):
            for f, file in enumerate(self.in_files):
                if reference_file == file:
                    reference_index = f
        else:
            SpecklepyTypeError('Reconstruction', 'reference_image', type(reference_file), 'int or str')

        return reference_index

    @property
    def reference_long_exp_file(self):
        return self.long_exp_files[self.reference_index]

    def select_box(self, radius):
        image = get_data(self.reference_long_exp_file, path=self.tmp_dir)
        logger.info(f"Select the source for the SSA reference aperture of radius {radius} pix!")
        _, selected = extract_sources(image=image, noise_threshold=3, fwhm=radius, select=True)
        if len(selected) == 0:
            sys.tracebacklimit = 0
            raise ValueError("No source selected to center the aperture on")
        pos = {'x': int(round(selected[0]['x'])), 'y': int(round(selected[0]['y']))}
        self.box = Box(indexes=[pos['x'] - radius, pos['x'] + radius + 1, pos['y'] - radius, pos['y'] + radius + 1])
        self.box.transpose()

    def create_long_exposures(self, integration_method=None, mask_hot_pixels=False, shifts=None):
        """Compute long exposures from the input data cubes.

        Arguments:
            integration_method (str, optional):
                Choose from 'collapse' or 'ssa' for either simply collapsing a cube along the time axis, or SSA'ing the
                individual frames.
            mask_hot_pixels (bool, optional):
                Let the SpeckleCube instance mask hot pixels.
            shifts (list, optional):
                List of shifts between the file and the reference file. This is used only in case of SSA
                reconstructions, where the reference box might require to be shifted.

        Returns:
            long_exposure_files (list):
                List of the file names of the long exposures.
        """

        # Update integration method, if provided
        if integration_method is not None:
            self.integration_method = integration_method

        # Initialize list of long exposure files
        long_exposure_files = []

        # Iterate over input data cubes
        for f, file in enumerate(self.in_files):

            # Read data from file
            speckle_cube = SpeckleCube(file_name=file, in_dir=self.in_dir, out_dir=self.tmp_dir,
                                       variance_extension=None)
            if speckle_cube.single_frame_mode:
                number_out_dir_levels = len(self.tmp_dir.split('/')) - 1
                long_exposure_files.append('../' * number_out_dir_levels + file)
                continue
            
            # Compute collapsed or SSA'ed images from the cube
            if self.integration_method == 'collapse':
                speckle_cube.collapse(mask=self.custom_mask)
            elif self.integration_method == 'ssa':
                if self.box is not None and shifts is not None:
                    box = self.box.shift(shift=shifts[f])
                else:
                    box = self.box
                box.crop_to_shape(speckle_cube.frame_shape)
                speckle_cube.ssa(box=box, mask_bad_pixels=mask_hot_pixels, mask=self.custom_mask)
            else:
                raise SpecklepyValueError('Reconstruction', 'integration_method', self.integration_method,
                                          expected="either 'collapse' or 'ssa'")

            # Mask out hot pixels
            # if mask_hot_pixels:
            #     speckle_cube.mask_hot_pixels()

            # Store data to a new Outfile instance
            logger.info(f"Saving temporary reconstruction of cube {file!r} to {speckle_cube.default_save_path()!r}")
            long_exposure_path = speckle_cube.write()

            # Add the recently created file to the list
            long_exposure_files.append(os.path.basename(long_exposure_path))

        return long_exposure_files

    def align_cubes(self, integration_method=None, alignment_mode='correlation', mask_hot_pixels=False):

        # Update integration method, if provided
        if integration_method is not None:
            self.integration_method = integration_method

        # Compute SSA reconstructions of cubes or collapse cubes for initial alignments
        if self.long_exp_files is None:
            self.long_exp_files = self.create_long_exposures(mask_hot_pixels=mask_hot_pixels)

        # Save time if only one cube is provided
        if self.single_cube_mode:
            self.shifts = [(0, 0)]
            self.pad_vectors = [((0, 0), (0, 0))]
            self.reference_pad_vector = ((0, 0), (0, 0))

        else:
            # Estimate relative shifts
            # self.shifts = alignment.estimate_shifts(files=self.long_exp_files, reference_file=self.reference_index,
            #                                         lazy_mode=True, return_image_shape=False, in_dir=self.tmp_dir,
            #                                         debug=self.debug)
            shift_estimator = alignment.ShiftEstimator()
            self.shifts = shift_estimator.estimate_shifts(file_names=self.long_exp_files, in_dir=self.tmp_dir,
                                                          reference_file_index=self.reference_index,
                                                          mode=alignment_mode, debug=self.debug)

            # Derive corresponding padding vectors
            self.pad_vectors, self.reference_pad_vector = \
                alignment.derive_pad_vectors(shifts=self.shifts, cube_mode=False, return_reference_image_pad_vector=True)

    def initialize_image(self):
        """Initialize the reconstruction image."""

        # Initialize image along the input frame shape
        image = np.zeros(self.frame_shape)

        # Crop or expand according to the reconstruction mode
        if self.mode == 'same':
            pass
        elif self.mode == 'full':
            image = np.pad(image, self.reference_pad_vector, mode='constant')
        elif self.mode == 'valid':
            # Estimate minimum overlap
            _shifts = np.array(self.shifts)
            _crop_by = np.max(_shifts, axis=0) - np.min(_shifts, axis=0)
            image = image[_crop_by[0]:, _crop_by[1]:]

        return image

    def coadd_long_exposures(self, integration_method=None, save=False):
        """Coadd the interim long exposures."""

        # Update integration method, if provided
        if integration_method is not None:
            self.integration_method = integration_method

        # Align cubes and cube relative shifts, if not available yet
        if self.shifts is None or self.pad_vectors is None or self.reference_pad_vector is None:
            self.align_cubes()

        # Initialize the image
        if self.image is None:
            self.image = self.initialize_image()

        # Iterate over long exposure images
        for index, file in enumerate(self.long_exp_files):

            # Read data
            tmp_image = get_data(file, path=self.tmp_dir)
            tmp_image_var = get_data(file, path=self.tmp_dir, extension=self.variance_extension,
                                     ignore_missing_extension=True)

            # Co-add reconstructions and var images
            self.image += alignment.pad_array(tmp_image, self.pad_vectors[index], mode=self.mode,
                                              reference_image_pad_vector=self.reference_pad_vector)
            if tmp_image_var is not None:
                if self.image_var is None:
                    self.image_var = self.initialize_image()
                self.image_var += alignment.pad_array(tmp_image_var, self.pad_vectors[index], mode=self.mode,
                                                      reference_image_pad_vector=self.reference_pad_vector)

        # Update out_file
        if save:
            self.save()

        return self.image, self.image_var

    def save(self):
        self.out_file.data = self.image
        if self.variance_extension is not None and self.image_var is not None:
            self.out_file.update_extension(ext_name=self.variance_extension, data=self.image_var)
