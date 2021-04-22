import glob
import os

from astropy.io import fits
from astropy.table import Table

from specklepy.logging import logger


def inspect(files, keywords, save=None, recursive=False, debug=False):

    # Find files
    if len(files) == 1:
        path = files[0]
        if '*' in path:
            files = glob.glob(path, recursive=recursive)
        else:
            files = glob.glob(os.path.join(path, '*fits'), recursive=recursive)

    # Terminal output
    if len(files):
        logger.info(f"Inspecting {len(files)} files for the FITS header keywords: {keywords}")
        files.sort()
    else:
        logger.error(f"Found no files in {files}!")
        raise RuntimeError(f"Found no files in {files}!")

    # Initialize output table
    out_table = Table(names=['FILE'] + keywords, dtype=[str] * (1 + len(keywords)))

    # Iterate through files and store header keywords
    for file in files:
        hdr = fits.getheader(filename=file)
        row = [file]
        for keyword in keywords:
            try:
                row.append(str(hdr[keyword]))
            except KeyError:
                row.append('--')
        out_table.add_row(row)

    # Store the table if requested
    if save is not None:
        out_table.write(save, format='ascii.fixed_width')

    # Present table
    print(out_table)
