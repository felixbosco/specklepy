def frame_shape(array):
    return array.shape[-2], array.shape[-1]


def frame_number(array):
    """Extract the number of frames in an array.

    Arguments:
        array (array-like):
            Array to extract the frame number from.

    Returns:
        frame_number (int):
            Number of frames in the array
    """
    if array.ndim == 2:
        number_frames = 1
    elif array.ndim == 3:
        number_frames = array.shape[0]
    else:
        raise ValueError(f"Data is {array.ndim}-dimensional but must be 2 or 3-dimensional!")

    return number_frames
