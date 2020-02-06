import numpy as np
import matplotlib.pyplot as plt

def get_sample_draw_vectors(nDraws, nFrames, first_uniform=False):
    """Generate draw vectors or weights for bootstrap resampling.

    Args:
        nDraws (int):
            Number of draws.
        nFrames (int):
            Number of frames to pick from.
        first_uniform (bool):
            Set True to receive the first draw as uniform wights/ all ones.

    Returns:
        sample_draw_vectors (numpy.ndarray):
            Array of shape (nDraws, nFrames) containing the counts of how often
            the particular frame at position (draw, frame) is picked.
    """

    shape = (nDraws, nFrames)
    draw_indizes = np.random.randint(nFrames, size=shape)

    sample_draw_vectors = np.zeros(shape=shape)
    for sample_draw_index, sample_draw_vector in enumerate(sample_draw_vectors):
        unique_indizes, counts = np.unique(draw_indizes[sample_draw_index], return_counts=True)
        sample_draw_vectors[sample_draw_index][unique_indizes] = counts

    if first_uniform:
        sample_draw_vectors[0] = 1

    return sample_draw_vectors


