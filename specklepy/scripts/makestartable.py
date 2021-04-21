import numpy as np

from astropy.units import Quantity

from specklepy.mock.startable import StarTable


def make_star_table(number_stars, iso_files, lf_files, lf_band, half_light_radius, population_weights=None,
                    out_file=None, table_format=None, overwrite=False):
    """Create a table of stars and write it to file

    Args:
        number_stars (int):
            Total number of stars to create.
        iso_files (list):
            List of isochrone files. Must at least contain one file name.
        lf_files (list):
            List of luminosity function files. Must contain the same number of file names as `iso_files`.
        lf_band (str):
            Name of the filter band that is represented in the luminosity function files. This is used for connecting
            the corresoonding tables. See documentation of `StarTable.generate_table` method for details.
        half_light_radius (int or float or str or Quantity):
            Half light radius of the generated cluster in linear units. Str-type values will be evaluated as a Quantity.
            Such input can be used to parse different units than parsecs.
        population_weights (list, optional):
            List of population weights, i.e. fraction of `number_stars` that is generated for the corresponding
            population, defined by the isochrone and luminosity function files. The weights are automatically normalized
            to unity. If not provided, all populations share the equal fraction of `number_stars`.
        out_file (str, optional):
            Name of a file, to which the star table shall be written.
        table_format (str, optional):
            Format string for writing the table to ASCII files. Falls back to `'ascii.fixed_width'` if not provided or
            to `'fits'` if `out_file` has the extension '.fits'.
        overwrite (bool, optional):
            Allow to overwrite existing table files.

    Returns:
        star_table (StarTable):
            StarTable instance that was created based on the provided isochrone and luminosity function files.
    """

    # Assert that the number of isochrone and luminosity function files are equal
    number_populations = len(iso_files)
    if len(lf_files) != number_populations:
        raise ValueError(f"Number of luminosity function files ({len(lf_files)}) does not match the number of "
                         f"isochrone files ({len(iso_files)})")

    # Derive population sizes
    if population_weights is None:
        population_weights = np.ones(number_populations)
    population_weights /= np.sum(population_weights)
    population_sizes_cum = np.cumsum(population_weights * number_stars).round().astype(int)
    population_sizes = np.append(population_sizes_cum[:1], np.diff(population_sizes_cum))

    # Initialize star table and generate the stellar populations
    star_table = StarTable()
    for population_size, iso_file, lf_file in zip(population_sizes, iso_files, lf_files):
        star_table.generate_table(population_size, iso_file, lf_file, lf_band, append=True)

    # Distribute stars in the field
    if isinstance(half_light_radius, str):
        half_light_radius = Quantity(half_light_radius)
    half_light_radius_unit = 'pc'
    if isinstance(half_light_radius, Quantity):
        half_light_radius_unit = half_light_radius.unit
        half_light_radius = half_light_radius.value
    star_table.add_spatial_columns(half_light_radius, half_light_radius_unit)

    # Write result to file
    if out_file is not None:
        if not out_file.endswith('.fits') and table_format is None:
            table_format = 'ascii.fixed_width'
        star_table.write(out_file, table_format=table_format, overwrite=overwrite)

    return star_table
