# wrappers
Various python/bash scripts wrapping around (mostly) astronomy-related software.

## sub_swarp.py
TODO: check fluxscale after `swarp` has done its thing. 

Use `swarp` (https://www.astromatic.net/software/swarp) to create a subimage from a FITS image. Useful particularly on large images where certain other methods fail (e.g. Montage's `mSubimage`) or simply take far longer than I would like. This can be used from the command line or by importing the script as module. Both methods use the `subimage` function. 
