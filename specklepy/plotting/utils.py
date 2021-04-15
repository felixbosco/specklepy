import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import os

from astropy.units import Quantity
from astropy.visualization import simple_norm

from specklepy.logging import logger
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.utils.box import Box


def imshow(image, title=None, norm=None, color_bar_label=None, save_to=None, maximize=False):
    """Shows a 2D image.

    Args:
        image (np.ndarray, ndim=2):
            Image to be plotted.
        title (str, optional):
            Plot title. Default is None.
        norm (str, optional):
            Can be set to 'log', for plotting in logarithmic scale. Default is
            None.
        color_bar_label (str, optional):
            Label of the color bar. Default is None.
        save_to (str, optional):
            Path to save the plot to. Default is None.
        maximize (bool, optional):
            Set true for showing the plot on full screen. Default is False.
    """

    if isinstance(image, np.ndarray):
        if image.ndim != 2:
            raise SpecklepyValueError('imshow()', 'image.ndim', image.ndim, '2')
        if isinstance(image, Quantity):
            unit = image.unit
            color_bar_label = "({})".format(unit)
            image = image.value
    else:
        raise SpecklepyTypeError('imshow()', 'image', type(image), 'np.ndarray')

    if norm == 'log':
        norm = clrs.LogNorm()
    else:
        norm = simple_norm(data=image, percent=99.)
    plt.figure()
    plt.imshow(image, norm=norm, origin='lower')
    plt.title(title)
    if maximize:
        maximize_plot()

    # Color bar
    color_bar = plt.colorbar(pad=0.0)
    if color_bar_label is not None:
        color_bar.set_label(color_bar_label)

    if save_to is not None:
        plt.savefig(save_to, dpi=300)

    plt.show()
    plt.close()


def plot3d(image):
    y, x = np.mgrid[-(image.shape[0] // 2): image.shape[0] // 2 + 1, -(image.shape[1] // 2): image.shape[1] // 2 + 1]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, image, cmap="viridis")
    plt.show()


def plot3d_aperture(image, pos, radius):
    box = Box.centered_at(pos[1], pos[0], radius)
    box.crop_to_shape(image.shape)
    aperture = box(image)
    plot3d(aperture)


def maximize_plot():
    mng = plt.get_current_fig_manager()
    backend = mpl.get_backend()
    if backend == 'Qt4Agg' or backend == 'Qt5Agg':
        mng.window.showMaximized()
    elif backend == 'wxAgg':
        mng.frame.Maximize(True)
    elif backend == 'TkAgg':
        try:
            mng.window.state('zoomed')
        except:
            logger.warning("Could not maximize plot")
    else:
        raise RuntimeWarning("Maximizing plot is not possible with matplotlib backend {}".format(backend))


def save_figure(file_name=None):
    """Save figure to a file, if a file name is provided.

    Args:
        file_name (str, optional):
            Name of the file to store the figure to. Nothing is done if not provided.
    """
    if file_name is not None:
        logger.info(f"Saving figure to {file_name}")

        # Identify requested file extension
        root, extension = os.path.splitext(file_name)
        extensions = ['.pdf', '.png', extension]

        # Save figure in multiple formats
        for ext in extensions:
            plot_file = root + ext
            try:
                plt.savefig(plot_file, bbox_inches='tight', pad_inches=0)
            except FileNotFoundError:
                path, file_name = os.path.split(plot_file)
                if not os.path.exists(path):
                    logger.info(f"Making dir {path}")
                    os.mkdir(path=path)
                plt.savefig(plot_file, bbox_inches='tight', pad_inches=0)


def desaturate_color(color, number_colors=1, saturation_values=None, saturation_min=0.1):
    """Desaturates a color and returns a list of desaturated colors.

    Args:
        color (str):
            Principle color to create desaturated representations of.
        number_colors (int, optional):
            Number of desaturated color representations.
        saturation_values (list, optional):
            List of saturation values.
        saturation_min (float, optional):
            Minimum value of saturation.

    Returns:
        colors (list):
            List of RGB `number_colors` desaturated representations of `color`.
    """

    # Input parameters
    if isinstance(color, str):
        rgb_color = clrs.to_rgb(color)
        hsv_color = clrs.rgb_to_hsv(rgb_color)
    elif isinstance(color, tuple):
        logger.info("Interpreting color tuple () as RGB values.")
        hsv_color = clrs.rgb_to_hsv(color)
    else:
        raise SpecklepyTypeError('desaturate_color()', 'color', type(color), 'str')

    if not isinstance(number_colors, int):
        raise SpecklepyTypeError('desaturate_color()', 'number_colors', type(number_colors), 'int')

    if not isinstance(saturation_min, float):
        raise SpecklepyTypeError('desaturate_color()', 'saturation_min', type(saturation_min), 'float')

    if saturation_values is None:
        saturation_values = np.linspace(hsv_color[1], saturation_min, num=number_colors)
    elif isinstance(saturation_values, list):
        pass  # list is correct, nothing to adapt
    elif isinstance(saturation_values, float):
        saturation_values = list(saturation_values)
    else:
        raise SpecklepyTypeError('desaturate_color()', 'saturation_values', type(saturation_values), 'list')

    # Create list of colors with varied saturation values
    colors = []
    for saturation_value in saturation_values:
        color = clrs.hsv_to_rgb((hsv_color[0], saturation_value, hsv_color[2]))
        colors.append(color)

    return colors
