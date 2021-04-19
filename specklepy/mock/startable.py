import numpy as np
from scipy.interpolate import interp1d

from astropy.table import Table


class StarTable(object):

    def __init__(self):
        self.table = None

    @staticmethod
    def sample_magnitudes_from_table(number_stars, table_file, table_format='ascii'):
        """Create samples of stellar magnitudes from a luminosity function table.

        Args:
            number_stars (int):
                Number of stars contained in the list.
            table_file (str):
                Name of the file that contains the tabulated luminosity function.
            table_format (str, optional):
                Format string for encoding the tabulated luminosity function.

        Returns:
            magnitudes (np.ndarray):
                Array of `number_stars` stellar magnitudes drawn from the luminosity function distribution.
        """
        table = Table.read(table_file, format=table_format)
        number_per_bin = np.power(10, table['col3'])
        to_magnitude = interp1d(x=number_per_bin, y=table['col2'], kind='cubic')
        draws = np.random.uniform(np.min(number_per_bin), np.max(number_per_bin), size=number_stars)
        return to_magnitude(draws)
