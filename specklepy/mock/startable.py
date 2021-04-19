import numpy as np
from scipy.interpolate import interp1d

from astropy.table import Table, QTable, Column

from specklepy.logging import logger


class StarTable(object):

    def __init__(self):
        self.table = None

    @property
    def length(self):
        return len(self.table)

    def __getitem__(self, item):
        return self.table.__getitem__(item)

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

        # Load luminosity function from table
        logger.info(f"Loading tabulated luminosity function from file {table_file!r}")
        table = Table.read(table_file, format=table_format)

        # Interpolate luminosity function
        logger.info(f"Interpolating luminosity function")
        number_per_bin = np.power(10, table['col3'])
        to_magnitude = interp1d(x=number_per_bin, y=table['col2'], kind='cubic')

        # Draw random numbers from uniform distribution and translate into magnitude space
        draws = np.random.uniform(np.min(number_per_bin), np.max(number_per_bin), size=number_stars)
        return to_magnitude(draws)

    def generate_table(self, number_stars, iso_table_file, lf_table_file, lf_band, table_format='ascii'):

        # Initialize aliases
        bad_columns = ['EEP', 'D51']
        units = {'EEP': None,
                 'M/Mo': 'solMass', 'LogTeff': 'Kelvin', 'LogG': None, 'LogL/Lo': 'solLum',
                 'U': 'mag', 'B': 'mag', 'V': 'mag', 'R': 'mag', 'I': 'mag',
                 'J': 'mag', 'H': 'mag', 'Ks': 'mag', 'Kp': 'mag',
                 'D51': None}

        # Initialize table
        self.table = QTable()

        # Draw samples of magnitudes in a given band
        magnitudes = self.sample_magnitudes_from_table(number_stars=number_stars, table_file=lf_table_file,
                                                       table_format=table_format)

        # Load data from iso-file
        logger.info(f"Loading tabulated isochrone data from file {iso_table_file!r}")
        iso_table = Table.read(iso_table_file, format=table_format)

        # Iterate through iso-table
        logger.info("Building table...")
        for column_name in iso_table.colnames:
            if column_name in bad_columns:
                continue
            elif column_name == lf_band:
                column = Column(data=magnitudes, name=column_name, unit='mag')
            else:
                # Interpolate tabulated functions
                to_quantity = interp1d(iso_table.columns[lf_band], iso_table.columns[column_name])

                # Insert transformed column of magnitudes into the current quantity
                column = Column(data=to_quantity(magnitudes), name=column_name, unit=units[column_name])
            logger.info(f"Inserting column for {column_name!r}")
            self.table.add_column(column)
