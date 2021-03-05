from datetime import datetime


def default_time_stamp():
    """Return a time stamp str of format 'YYYYMMDD_HHMMSS'."""
    return datetime.now().strftime('%Y%m%d_%H%M%S')
