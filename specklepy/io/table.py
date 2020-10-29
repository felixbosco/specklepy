from astropy.io.registry import IORegistryError
from astropy.table import Table


def read_table(file_name, format=None, error=True):

    # Initialize output table
    table = None

    # Read table from file
    if format is not None:
        table = Table.read(file_name, format=format)
    else:
        # Try table default table formats
        for _format in ['ascii.fixed_width', 'ascii', None]:
            try:
                table = Table.read(file_name, format=_format)
            except IORegistryError:
                pass

    if error and table is None:
        raise RuntimeError(f"Reading the table from {file_name} was not successful!")

    return table
