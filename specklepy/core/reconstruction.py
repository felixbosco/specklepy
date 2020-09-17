import numpy as np
import os

from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.logging import logger


class Reconstruction(object):

    """Base class for image reconstructions.

    This class stores the principle parameters of the reconstructed image. In the future, the functions `ssa` and
    `holography` may become child classes of this one.
    """

    supported_modes = ['full', 'same', 'valid']

    def __init__(self, in_files, mode='full', reference_image=None, out_file=None, in_dir=None):
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
            reference_image (int or str, optional):
                The index in the `in_files` list or the name of the image serving as reference in 'same' mode.
            out_file (str, optional):
                Name of an output file to store the reconstructed image in.
            in_dir (str, optional):
                Path to the `in_files`.
        """

        # Check input parameter types
        if not isinstance(in_files, list):
            raise SpecklepyTypeError('Reconstruction', 'in_files', type(in_files), 'list')
        if not isinstance(mode, str):
            raise SpecklepyTypeError('Reconstruction', 'mode', type(mode), 'str')
        if out_file is not None and not isinstance(out_file, str):
            raise SpecklepyTypeError('Reconstruction', 'outfile', type(out_file), 'str')
        if reference_image is not None and not isinstance(reference_image, (int, str)):
            raise SpecklepyTypeError('Reconstruction', 'reference_image', type(reference_image), 'int or str')

        # Check input parameter values
        if mode not in self.supported_modes:
            raise SpecklepyValueError('Reconstruction', 'mode', mode, f"in {self.supported_modes}")
        if mode is 'same' and reference_image is None:
            logger.warning("Reconstruction requires a reference image if run in `same` mode! The first image will "
                           "be used as reference field.")
            reference_image = 0
        if in_dir is None:
            in_dir = ''

        # Store input data
        self.in_files = in_files
        self.mode = mode
        self.outfile = out_file
        self.reference_image = reference_image

        # Derive shape of individual input frames
        single_cube_mode = len(self.in_files) == 1
        example_frame = fits.getdata(os.path.join(in_dir, self.in_files[0]))
        if example_frame.ndim == 3:
            example_frame = example_frame[0]
        frame_shape = example_frame[0]

        # Initialize image
        if single_cube_mode:
            self.image = np.zeros(frame_shape)
            self.shifts = (0, 0)
        else:
            # Estimate relative shifts
            self.shifts = None

            # Derive corresponding padding vectors
            pass

            # Derive corresponding image sizes
            if self.mode == 'same':
                self.image = np.zeros(frame_shape)
            elif self.mode == 'full':
                pass
            elif self.mode == 'valid':
                pass

        # Initialize out file
        pass
