import numpy as np
import os
from astropy.io import fits
from datetime import datetime

from specklepy.exceptions import SpecklepyTypeError
from specklepy.io.outfile import Outfile
from specklepy.logging import logger


class ReductionFile(Outfile):

    def __init__(self, file, data=None, prefix=None, path=None, reduction=None, last_reduction=None, header_card_prefix="HIERARCH SPECKLEPY"):
        """Class that carries the link to a file for data reduction products.

        Args:
            file (str):
                Input file that will be serve as a source for the header and
                raw data.
            data (np.ndarray, optional):
                Data that will be stored as the new PrimaryHDU.
            prefix (str, optional):
                Prefix that will combined with the 'file' argument to the name
                of the new file.
            path (str, optional):
                Target path under which the new file will be stored.
            reduction (str, optional):
                Current reduction step, will be stored to the file header.
            last_reduction (str, optional):
                Last reduction step, will be used as a name for the new fits extension.
            header_card_prefix (str, optional):
                Default prefix for header cards.
        """

        # Check input parameters
        if isinstance(file, str):
            self.parent_file = file
        else:
            raise SpecklepyTypeError('ReductionFile', argname='file', argtype=type(file), expected='str')

        if data is None or isinstance(data, np.ndarray):
            self._data = data
        else:
            raise SpecklepyTypeError('ReductionFile', argname='data', argtype=type(data), expected='np.ndarray')

        if prefix is None:
            self.prefix = ""
        elif isinstance(prefix, str):
            self.prefix = prefix
        else:
            raise SpecklepyTypeError('ReductionFile', argname='prefix', argtype=type(prefix), expected='str')

        if path is None or isinstance(path, str):
            self.path = path
        else:
            raise SpecklepyTypeError('ReductionFile', argname='path', argtype=type(path), expected='str')

        if reduction is None or isinstance(reduction, str):
            self.reduction = reduction
        else:
            raise SpecklepyTypeError('ReductionFile', argname='reduction', argtype=type(reduction), expected='str')

        if last_reduction is None:
            self.last_reduction = "RAW"
        elif isinstance(last_reduction, str):
            self.last_reduction = last_reduction
        else:
            raise SpecklepyTypeError('ReductionFile', argname='last_reduction', argtype=type(last_reduction), expected='str')

        if header_card_prefix is None or isinstance(header_card_prefix, str):
            self.header_card_prefix = header_card_prefix
        else:
            raise SpecklepyTypeError('ReductionFile', argname='header_card_prefix', argtype=type(header_card_prefix), expected='str')


        # Create file name
        self.filename = self.prefix + os.path.basename(self.parent_file) # Make sure to get rid of the path


        # Read header information and data from parent file
        with fits.open(self.parent_file) as hdulist:
            # Copy parent file data into extensions
            extensions = []
            for hdu in hdulist:
                ext = {'name': hdu.name,
                       'data': hdu.data,
                       'header': hdu.header}
                extensions.append(ext)

            # Construct primary HDU
            header = hdulist[0].header

            if 'PIPELINE' not in header.keys():
                # Parent file is not a ReductionFile
                header.set('PIPELINE', 'SPECKLEPY')
            header.set(reduction, str(datetime.now())) # Store the current reduction step

            if self._data is None:
                primary_data = hdulist[0].data
            else:
                primary_data = self._data
            del self._data

        # parent_data, header = fits.getdata(self.parent_file, header=True)

        # if self._data is None:
        #     primary_hdu = fits.PrimaryHDU(data=parent_data, header=header)
        # else:
        #     primary_hdu = fits.PrimaryHDU(data=self._data, header=header)
        # hdulist = fits.HDUList(hdus=[primary_hdu])
        # parent_image_hdu = fits.ImageHDU(data=parent_data, name=self.last_reduction)
        # hdulist.append(parent_image_hdu)
        # hdulist.writeto(os.path.join(self.path, self.filename))


        # todo check whether this is working in flat!
        super().__init__(filename=self.filename,
                         path=self.path,
                         data=primary_data,
                         extensions=extensions,
                         header=header,
                         cards=None,
                         timestamp=False,
                         header_card_prefix=header_card_prefix)