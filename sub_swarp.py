#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import uuid
import shutil
import subprocess32 as subprocess
from astropy.io import fits


def config_file_formatter(outname, ra, dec, pscale, imsize, projection="SIN", 
                          memmax=4095):
    """Create configuration file for swarp to create a subimage.

    There are additional options for swarp, but most are for co-addition.
    """


    with open(outname.replace(".fits", "")+".swarp", "w+") as f:
        f.write("""
# Default configuration file for SWarp 2.38.0
# EB 2018-06-02
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME          {0}.fits        # Output filename
WEIGHTOUT_NAME       {0}.weight.fits   # Output weight-map filename
 
HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers
  
#-------------------------------- Astrometry ----------------------------------
 
CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        {1}             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            MANUAL          # MANUAL, ALL or MOST
CENTER               {2}, {3}          # Coordinates of the image center
PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            {4}             # Pixel scale (arcsec)
IMAGE_SIZE             {5}             # Image size (0 = AUTOMATIC) (pixels)
 
#-------------------------------- Resampling ----------------------------------
 
RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images
 
RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           4               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)
 
#------------------------------ Memory management -----------------------------
 
VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               {6}             # Maximum amount of virtual memory (MB)
MEM_MAX                {6}             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        {6}             # RAM dedicated to co-addition(MB)
 
""".format(outname.replace(".fits", ""), projection, ra, dec, pscale*3600., 
           int(imsize), memmax))

    return outname.replace(".fits", "")+".swarp"



def bpp(major, minor, cd1, cd2):
    """Conversion from Jy/beam to Jy/pixel."""

    return np.pi*(major*minor) / (cd1*cd2 * 4.*np.log(2.))



def subimage(image, imsize, outname=None, ra=None, dec=None, units="deg", 
             projection="native"):
    """Create subimage with swarp.

    First generate a configuration file, then run swarp.
    
    Parameters
    ----------
    image : str
        FITS image to create subimage from.
    imsize : float
        Size of image in either degrees [`units="deg"`] or 
        pixels [`units="pixels"`].
    outname : str, optional
        Output file name [image_sub.fits].
    ra : float, optional
        RA centre of subimage in degrees. [Default is CRVAL1 from `image`].
    dec : float, optional
        DEC centre of subimage in degrees. [Default is CRVAL2 from `image`].
    units : str, ["degrees" or "pixels"]
        Units of `imsize`. ["deg"]. 
    projection : str, optional
        Projection for output image. See swarp documentation for options. 
        "native" will give the same projection as input image. ["native"].
    
    """

    hdr = fits.getheader(image)
    if "CDELT1" in hdr.keys():
        pscale = abs(hdr["CDELT1"])
        cd2 = abs(hdr["CDELT2"])
    elif "CD1_1" in hdr.keys():
        pscale = abs(hdr["CD1_1"])
        cd2 = abs(hdr["CD2_2"])
    else:
        raise ValueError("No CDELT1 or CD1_1 found.")

    if ra is None:
        ra = hdr["CRVAL1"]
    if dec is None:
        dec = hdr["CRVAL2"]

    if "d" in units.lower():
        imsize /= pscale

    if "native" in projection.lower():
        projection = hdr["CTYPE1"].split("-")[-1]

    # Swarp will strip beam information:
    beam = None
    if "BMAJ" in hdr.keys():
        beam = (hdr["BMAJ"], hdr["BMIN"], hdr["BPA"])
    if beam is not None:
        # We need to carefully preserve "Jy/beam" units by sripping them,
        # then adding them back after making the subimage.
        # This should only be necessary if the pixel scale changes, which it 
        # shouldn't if only making a subimage like this.
        # HOWEVER, I was caught out when resampling to different pixel scales,
        # so this is here in case the resampling causes issues.
        b = bpp(major=beam[0],
                minor=beam[1],
                cd1=pscale,
                cd2=cd2)

        # First generate a temporary file name:
        tempf = str(uuid.uuid4())+".fits"
        shutil.copy2(image, tempf)
        # Then convert from Jy/beam to Jy/pixel:
        with fits.open(tmpf, mode="update") as tf:
            tf[0].data /= b
            tf.flush()

        sub_image = tempf
    else:
        sub_image = image

    # Frequency data might not be preserved either:
    freq = None
    if "FREQ" in hdr.keys():
        freq = hdr["FREQ"]
    elif "CRVAL3" in hdr.keys():
        freq = hdr["CRVAL3"]

    if outname is None:
        outname = image.replace(".fits", "_sub.fits")

    # Generate config. file:
    cfg = config_file_formatter(outname, ra, dec, pscale, imsize, projection)
    s = "swarp {0} -c {1}".format(sub_image, cfg)
    subprocess.call(s, shell=True)


    # Add beam and/or frequency key back:
    with fits.open(outname, mode="update") as f:
        if beam is not None:
            f[0].header["BMAJ"] = beam[0]
            f[0].header["BMIN"] = beam[1]
            f[0].header["BPA"] = beam[2]
            b = bpp(major=beam[0],
                    minor=beam[1],
                    cd1=abs(f[0].header["CDELT1"]),
                    cd2=abs(f[0].header["CDELT2"]))
            f[0].data *= b

        if freq is not None:
            f[0].header["FREQ"] = freq

        f.flush()

    if beam is not None:
        os.remove(tempf)



# ---------------------------------------------------------------------------- #
if __name__ == "__main__":

    import argparse
    ps = argparse.ArgumentParser(description="Create a subimage of a FITS file " 
                                             "using SWARP.",
                                 epilog="If multiple images are given, either "
                                        "multiple outnames should also be given "
                                        "otherwise outnames will follow the "
                                        "default format.")
    ps.add_argument("-i", "--image", dest="image", nargs="+",
                    help="FITS image to create subimage from.", required=True)
    ps.add_argument("-s", "--imsize", dest="imsize", help="Size of image in "
                    "either degrees or pixels (with --pixels).",
                    type=float, required=True)
    ps.add_argument("-o", "--outname", dest="outname", nargs="+", 
                    help="Output file name [image_sub.fits]", default=[None])
    ps.add_argument("-r", "--ra", dest="ra", help="RA centre of subimage in "
                    "degrees. [Default is CRVAL1 from `image`].", type=float)
    ps.add_argument("-d", "--dec", dest="dec", help="DEC centre of subimage in "
                    "degrees. [Default is CRVAL2 from `image`].", type=float)
    ps.add_argument("--pixels", dest="unit", const="pixels", action="store_const",
                    help="Switch corresponding to `imsize` being specified "
                    " in units of image pixels.", default="deg")
    ps.add_argument("-p", "--projection", dest="projection", help="Projection "
                    "for output image. See swarp documentation for options. "
                    "'native' will give the same projection as input image. "
                    "['native'].", default="native")
    args = ps.parse_args()

    if len(args.outname) != len(args.image):
        args.outname = [None]*len(args.image)
    
    for i, image in enumerate(args.image):
        subimage(image, args.imsize, args.outname[i], args.ra, args.dec,
                 args.unit, args.projection) 
        print("Subimage of {0} made.".format(image))
    sys.exit(0)




