from IPython import embed
import numpy as np
import os

from astropy.io import fits

from specklepy.core import alignment
from specklepy.core.specklecube import SpeckleCube
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.reconstructionfile import ReconstructionFile
from specklepy.logging import logger
from specklepy.utils.box import Box


class Reconstruction(object):

    """Base class for image reconstructions.

    This class stores the principle parameters of the reconstructed image. In the future, the functions `ssa` and
    `holography` may become child classes of this one.
    """

    supported_modes = ['full', 'same', 'valid']

    def __init__(self, in_files, mode='same', reference_file=None, out_file=None, in_dir=None, tmp_dir=None,
                 integration_method='collapse', var_ext=None, box_indexes=None, debug=False):
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
        self.var_ext = var_ext  # if var_ext is not None else 'VAR'
        self.box = Box(box_indexes) if box_indexes is not None else None
        self.debug = debug

        # Derive shape of individual input frames
        example_frame = fits.getdata(os.path.join(in_dir, self.in_files[0]))
        if example_frame.ndim == 3:
            example_frame = example_frame[0]
        self.frame_shape = example_frame.shape

        # Initialize secondary attributes
        self.image = None
        self.image_var = None
        self.long_exp_files = None
        self.shifts = None
        self.pad_vectors = None
        self.reference_pad_vector = None

        # Initialize image
        if self.single_cube_mode:
            self.image = np.zeros(self.frame_shape)
            self.shifts = (0, 0)
        else:
            # Align cubes
            self.align_cubes()

            # Derive corresponding image sizes
            self.image = self.initialize_image()

        # Initialize the variance map
        self.image_var = np.zeros(self.image.shape) if self.var_ext is not None else None

        # Initialize output file and create an extension for the variance
        self.out_file = ReconstructionFile(files=self.in_files, filename=self.out_file, shape=self.image.shape,
                                           in_dir=in_dir, cards={"RECONSTRUCTION": "SSA"})
        if self.image_var is not None:
            self.out_file.new_extension(name=self.var_ext, data=self.image_var)

    @property
    def single_cube_mode(self):
        return len(self.in_files) == 1

    @property
    def reference_file(self):
        return self.in_files[self.reference_index]

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

    def create_long_exposures(self, integration_method=None):
        """Compute long exposures from the input data cubes."""

        # Update integration method, if provided
        if integration_method is not None:
            self.integration_method = integration_method

        # Initialize list of long exposure files
        long_exposure_files = []

        # Iterate over input data cubes
        for file in self.in_files:

            # Read data from file
            speckle_cube = SpeckleCube(file_name=file, in_dir=self.in_dir, out_dir=self.tmp_dir,
                                       variance_extension=None)

            # Compute collapsed or SSA'ed images from the cube
            if self.integration_method == 'collapse':
                speckle_cube.collapse()
            elif self.integration_method == 'ssa':
                speckle_cube.ssa(box=self.box)
            else:
                raise SpecklepyValueError('Reconstruction', 'integration_method', self.integration_method,
                                          expected="either 'collapse' or 'ssa'")

            # Store data to a new Outfile instance
            logger.info(f"Saving temporary reconstruction of cube {file} to {speckle_cube.default_save_path()}")
            long_exposure_path = speckle_cube.store()

            # Add the recently created file to the list
            long_exposure_files.append(os.path.basename(long_exposure_path))

        return long_exposure_files

    def align_cubes(self, integration_method=None):

        # Update integration method, if provided
        if integration_method is not None:
            self.integration_method = integration_method

        # Compute SSA reconstructions of cubes or collapse cubes for initial alignments
        self.long_exp_files = self.create_long_exposures()

        # Estimate relative shifts
        self.shifts = alignment.get_shifts(files=self.long_exp_files, reference_file=self.reference_index,
                                           lazy_mode=True, return_image_shape=False, in_dir=self.tmp_dir,
                                           debug=self.debug)

        # Derive corresponding padding vectors
        self.pad_vectors, self.reference_pad_vector = \
            alignment.get_pad_vectors(shifts=self.shifts, cube_mode=False, return_reference_image_pad_vector=True)

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

    def coadd_long_exposures(self, integration_method='ssa', save=False):
        """Coadd the interim long exposures."""

        # Update integration method, if provided
        if integration_method is not None:
            self.integration_method = integration_method

        # Align cubes and cube relative shifts, if not available yet
        if self.shifts is None or self.pad_vectors is None or self.reference_pad_vector is None:
            self.align_cubes()

        # Iterate over long exposure images
        for index, file in enumerate(self.long_exp_files):

            # Read data
            with fits.open(os.path.join(self.tmp_dir, file)) as hdu_list:
                tmp_image = hdu_list[0].data
                if self.var_ext is not None and self.var_ext in hdu_list:
                    tmp_image_var = hdu_list[self.var_ext].data
                else:
                    tmp_image_var = None

            # Co-add reconstructions and var images
            self.image += alignment.pad_array(tmp_image, self.pad_vectors[index], mode=self.mode,
                                              reference_image_pad_vector=self.reference_pad_vector)
            if tmp_image_var is not None:
                self.image_var += alignment.pad_array(tmp_image_var, self.pad_vectors[index], mode=self.mode,
                                                      reference_image_pad_vector=self.reference_pad_vector)

        # Update out_file
        if save:
            self.save()

        return self.image, self.image_var

    def save(self):
        self.out_file.data = self.image
        if self.var_ext is not None and self.image_var is not None:
            self.out_file.update_extension(ext_name=self.var_ext, data=self.image_var)
