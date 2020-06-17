class Reconstruction(object):

    """Base class for image reconstructions.

    This class stores the principle parameters of the reconstructed image. In the future, the functions `ssa` and
    `holography` may become child classes of this one.
    """

    def __init__(self, in_files, mode='full', outfile=None, reference_image=None):
        """Create a Reconstruction instance.

        Args:
            in_files (list):
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
            outfile (str, optional):
                Name of an output file to store the reconstructed image in.
            reference_image (int or str, optional):
                The index in the `in_files` list or the name of the image serving as reference in 'same' mode.
        """

        # Store input data
        self.in_files = in_files
        self.mode = mode
        self.outfile = outfile
        self.reference_image = reference_image

        # Initialize image
        pass