"""
montage_reproj
==============

Re-project a set of FITS images
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np
import array, time, sets

import matplotlib.pyplot as plt
import glob

from astropy.table import Table


def main(fitsfile=None, imgdir=None, silent=True, verbose=False,
         catalog='SDSS'):
    '''
    Main function for montage_reproj

    Parameters
    ----------
    fitsfile : string
      Filename of FITS file containing images in different FITS extensions
      Default: None. Either specify this or imgdir, but not both

    imgdir: string
      Full path to where images are located. Default: None.
      Either sspecify this or fitsfile, but not both

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    catalog : string
      Which survey (e.g., SDSS, 2MASS) the finding chart image is from.
      Default: 'SDSS'
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 2 January 2017
    '''

    if silent == False:
        print '### Begin montage_reproj.main | '+systime()

    if fitsfile != None:
        imgdir = '/var/tmp/montage/'
        if not exists(imgdir): os.mkdir(imgdir)

        hdu = fits.open(fitsfile)
        n_images = len(hdu)

        if n_images == 1:
            print '## Only one image is provided'
        else:
            for ii in range(n_images):
                hdu[ii].writeto(imgdir+'im'+str(ii)+'.fits', clobber=True)

    imgtable = imgdir+'im.tbl'
    montage.mImgtbl(imgdir, imgtable)

    if silent == False:
        print '### End montage_reproj.main | '+systime()
#enddef
