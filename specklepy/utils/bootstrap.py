import numpy as np


def random_draw_vectors(number_draws, number_frames, first_uniform=False):
    """Generate draw or weight vectors for bootstrap resampling.

    Args:
        number_draws (int):
            Number of draws.
        number_frames (int):
            Number of frames to pick from.
        first_uniform (bool):
            Set True to receive the first draw as uniform wights/ all ones.

    Returns:
        sample_draw_vectors (numpy.ndarray):
            Array of shape (nDraws, nFrames) containing the counts of how often
            the particular frame at position (draw, frame) is picked.
    """

    # Draw the frame indexes randomly
    shape = (number_draws, number_frames)
    draw_indizes = np.random.randint(number_frames, size=shape)

    # Combine into array of vectors with weights ("counts")
    sample_draw_vectors = np.zeros(shape=shape)
    for sample_draw_index, sample_draw_vector in enumerate(sample_draw_vectors):
        unique_indizes, counts = np.unique(draw_indizes[sample_draw_index], return_counts=True)
        sample_draw_vectors[sample_draw_index][unique_indizes] = counts

    # Replace first draw vector with a homogeneous one
    if first_uniform:
        sample_draw_vectors[0] = 1

    return sample_draw_vectors
