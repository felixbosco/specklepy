class Segmentation(object):

    """Segmentation of fields of view.

    This class handles the segmentation of fields of view in (nx, ny) segments.

    Attributes:
        nx (int):
            Number of segments along x-axis (0-index).
        ny (int):
            Number of segments along y-axis (1-index).

    """

    def __init__(self, nx, ny, image_shape):
        """Create a Segmentation instance.

        Args:
            nx (int):
                Number of segments along x-axis (0-index).
            ny (int):
                Number of segments along y-axis (1-index).
            image_shape (tuple):
                Shape of the image to be segmented.
        """

        # Store parameters
        self.nx = nx
        self.ny = ny
        if len(image_shape) == 2:
            self.image_shape = image_shape
        elif len(image_shape) == 3:
            # First axis is assumed to be time axis.
            self.image_shape = tuple([image_shape[1], image_shape[2]])