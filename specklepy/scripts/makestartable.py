import numpy as np

from specklepy.mock.startable import StarTable


def make_star_table(number_stars, iso_files, lf_files, lf_band, out_file, population_weights=None):

    # Assert that the number of isochrone and luminosity function files are equal
    number_populations = len(iso_files)
    if len(lf_files) != number_populations:
        raise ValueError(f"Number of luminosity function files ({len(lf_files)}) does not match the number of "
                         f"isochrone files ({len(iso_files)})")
    if population_weights is None:
        population_weights = np.ones(number_populations)

    star_table = StarTable()
    star_table.generate_table(number_stars, iso_file, lf_file, lf_band)
    star_table.add_spatial_columns(27, 'pc')
    star_table.write(out_file)
    star_table.write(out_file.replace('.dat', '.fits'), table_format='fits')
