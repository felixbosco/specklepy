# Specklepy - A versatile tool for analyzing astronomical short-exposure ("speckle") data

## To Do:
* Make logger write a copy to a .log file
* Implement RAM save mode, which saves the Fourier transformed images and PSFs to files and then just sums up within a tmp file
* Create general script with signature "specklepy holography -f ..."
* Tune star finder routines
* [low priority] implement 'valid' mode, that yields a reconstruction only in the intersection of all file field of views
* Implement data reduction codes
* Implement frame generation via importing vegapy
  * Implement dithering
* For holography() implement secondary source subtraction
