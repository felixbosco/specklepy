import numpy as np


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

        # Compute widths of the segments
        self.dx = self.image_shape[0] // self.nx
        self.dy = self.image_shape[1] // self.ny

        # Derive a list of field segments
        self.segments = []
        for ix in range(self.nx):
            for iy in range(self.ny):
                xmin = ix * self.image_shape[0] // self.nx
                xmax = (ix+1) * self.image_shape[0] // self.nx - 1
                ymin = iy * self.image_shape[1] // self.ny
                ymax = (iy+1) * self.image_shape[1] // self.ny - 1

                self.segments.append(Segment(xmin, xmax, ymin, ymax))

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        try:
            tmp = self.segments[self.index]
        except IndexError:
            raise StopIteration
        self.index += 1
        return tmp

    def __getitem__(self, item):
        return self.segments[item]

    def all_covered(self, positions):
        """Are all segments in this Segmentation covered by at least one member of a position list?

        Args:
            positions (list of tuple):
                List of 2-tuples of positions.

        Returns:
            all_covered (bool):
                True, if each segment is covered by at least one element of positions.
        """

        covered = np.zeros(len(self.segments), dtype=bool)

        for pos in positions:
            for s, seg in enumerate(self.segments):
                covered[s] |= pos in seg

        return covered.all()


class Segment(object):

    """Class for storing the properties of individual segments in an image segmentation.

    Attributes:
        xmin (int):
            Minimum index of the segment within the parental field, along x-axis.
        xmax (int):
            Maximum index of the segment within the parental field, along x-axis.
        ymin (int):
            Minimum index of the segment within the parental field, along y-axis.
        ymax (int):
            Maximum index of the segment within the parental field, along y-axis.
    """

    def __init__(self, xmin, xmax, ymin, ymax):
        """Create a Segment object.

        Args:
            xmin (int):
                Minimum index of the segment within the parental field, along x-axis.
            xmax (int):
                Maximum index of the segment within the parental field, along x-axis.
            ymin (int):
                Minimum index of the segment within the parental field, along y-axis.
            ymax (int):
                Maximum index of the segment within the parental field, along y-axis.
        """

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def __str__(self):
        return f"Segment: [{self.xmin}:{self.xmax}, {self.ymin}:{self.ymax}]"

    def __call__(self, array):
        """Return the segment of an image array.

        Args:
            array (array_like):
                Image array from which the segment shall be extracted.
        """
        return array[self.xmin:self.xmax, self.ymin:self.ymax]

    def __contains__(self, pos):
        """Return bool of whether a position tuple lies within the segment boundaries.

        Args:
            pos (tuple):
                Position tuple of (x, y) coordinates.

        Returns:
            in (bool):
                True if pos lies within segment limits.
        """
        return self.xmin <= pos[0] <= self.xmax and self.ymin <= pos[1] <= self.ymax
