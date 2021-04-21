import numpy as np

from astropy.units import Quantity

from specklepy.mock.startable import StarTable


def make_star_table(number_stars, iso_files, lf_files, lf_band, out_file, half_light_radius, population_weights=None):

    # Assert that the number of isochrone and luminosity function files are equal
    number_populations = len(iso_files)
    if len(lf_files) != number_populations:
        raise ValueError(f"Number of luminosity function files ({len(lf_files)}) does not match the number of "
                         f"isochrone files ({len(iso_files)})")

    # Derive population sizes
    if population_weights is None:
        population_weights = np.ones(number_populations)
    population_weights /= np.sum(population_weights)
    population_sizes = np.empty(number_populations, dtype=int)
    for w, weight in enumerate(population_weights):
        population_sizes[w] = int(weight * number_stars)
    while np.sum(population_sizes) < number_stars:
        population_sizes[np.argmax(population_weights)] += 1

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
    star_table.write(out_file)
