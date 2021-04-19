import numpy as np
from scipy.interpolate import interp1d

from astropy.table import Table, QTable, Column
from astropy.units import Unit, UnitConversionError, Quantity

from specklepy.logging import logger


class StarTable(object):

    def __init__(self):
        self.table = None

    @classmethod
    def from_file(cls, file_name, table_format=None):
        """Initialize a StarTable instance from a table file.

        Args:
            file_name (str):
                Name of or path to the file containing the steller parameters.
            table_format (str, optional):
                 Format string for encoding the tabulated luminosity function.

        Returns:
            obj (StarTable):
                StarTable object as initialized from the table file.
        """
        obj = cls.__init__()
        obj.table = Table.read(file_name, format=table_format)
        return obj

    def __len__(self):
        return len(self.table)

    def __getitem__(self, item):
        return self.table.__getitem__(item)

    def __repr__(self):
        return self.table.__repr__()

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
        """Generate a complete table for all the quantities represented in a given isochrone file.

        Args:
            number_stars (int):
                Number of stars to generate.
            iso_table_file (str):
                Name of or path to the file containing tabulated isochrones.
            lf_table_file (str):
                Name of or path to the file containing tabulated luminosity functions.
            lf_band (str):
                Name of the filter band, which is represented in the luminosity function file. This is used for linking
                these magnitudes to the other quantities from the isochrone table.
            table_format (str, optional):
                 Format string for encoding the tabulated luminosity function.
        """

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

    def add_spatial_columns(self, half_light_radius, half_light_radius_unit='pc'):
        """Add three spatial coordinates to the stars that are normally distributed about the cluster center.

        Args:
            half_light_radius (float-like):
                Value of the star clusters half light radius. The units can be set via the argument.
            half_light_radius_unit (str, optional):
                Str-type representation of a unit.
        """
        
        # Derive standard deviation in radial direction from half-light radius
        radial_scale = half_light_radius * 2 / 2.35

        # Iterate through spatial axes
        for column_name in ['x', 'y', 'z']:
            pos = np.random.normal(loc=0.0, scale=radial_scale, size=len(self))
            column = Column(name=column_name, data=pos, unit=half_light_radius_unit)
            logger.info(f"Inserting column for {column_name!r}")
            self.table.add_column(column)

    def transfer_column_to_angles(self, key, distance, new_name=None, default_linear_unit='pc', output_unit='arcsec'):

        # Get column from table
        column_data = self.table[key]

        # Assign default unit
        if column_data.unit is None:
            logger.warning(f"Assuming that x-axis data are provided in units of parsecs!")
            column_data.unit = default_linear_unit

        # Try to convert units
        try:
            column_data = np.arctan(np.divide(column_data, distance))
        except UnitConversionError:
            logger.warning(f"Unable to convert from units of {column_data.unit!r} to {default_linear_unit!r}!")

        # Assign default to `new_name`
        if new_name is None:
            new_name = key

        return Column(name=new_name, data=column_data.to(output_unit))

    @staticmethod
    def distance_modulus(distance):
        return 5 * np.log10(np.divide(distance, Quantity('10 pc'))) * Unit('mag')

    def apply_distance_modulus(self, key, distance, new_name=None, data_unit='mag', output_unit='mag'):

        # Get column from table
        column_data = self.table[key]

        # Transform from any units to magnitudes
        if data_unit != 'mag':
            raise NotImplementedError(f"Handling of units other 'mag' is not implemented yet! (Received {data_unit})")

        # Add distance modulus
        try:
            column_data += self.distance_modulus(distance=distance)
        except UnitConversionError:
            column_data = column_data * Unit('mag') + self.distance_modulus(distance=distance)

        # Assign default to `new_name`
        if new_name is None:
            new_name = key

        return Column(name=new_name, data=column_data.to(output_unit))

    def observe_at_distance(self, distance, filter_band, x_key='x', y_key='y', filter_band_unit='mag',
                            distance_unit='pc'):

        # Add default unit to distance parameters
        if not hasattr(distance, 'unit'):
            logger.warning(f"Assuming default unit of {distance_unit!r} for distance without unit!")
            distance = distance * Unit(distance_unit)

        # Initialize output table
        out_table = Table()

        # Scale positional coordinates
        out_table.add_column(self.transfer_column_to_angles(key=x_key, distance=distance))
        out_table.add_column(self.transfer_column_to_angles(key=y_key, distance=distance))

        # Scale photometry
        out_table.add_column(self.apply_distance_modulus(key=filter_band, distance=distance,
                                                         data_unit=filter_band_unit))

        return out_table

    def write(self, file_name, table_format='ascii.fixed_width', overwrite=True):
        """Save the table to an ASCII or FITS file.

        Args:
            file_name (str):
                Name or path of the file to save the table to
            table_format (str, optional):
                Format of the
            overwrite (bool, optional):
                Allow to overwrite existing files. Default is True.
        """
        logger.info(f"Writing table data to file {file_name!r}")
        self.table.write(file_name, format=table_format, overwrite=overwrite)
