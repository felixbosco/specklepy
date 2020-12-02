import numpy as np

from tqdm import trange

from astropy import stats


def sigma_clipped_stats(cube, **kwargs):
    number_rows = cube.shape[1]
    frame_shape = cube.shape[1:]
    mean = np.empty(frame_shape)
    median = np.empty(frame_shape)
    std = np.empty(frame_shape)

    for row_index in trange(number_rows[:2]):
        row_stats = stats.sigma_clipped_stats(cube[:, row_index, :], axis=0, **kwargs)
        mean[row_index], median[row_index], std[row_index] = row_stats

    return mean, median, std
