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
import glob

import astropy.units as u
from astropy.wcs import WCS
from astropy import wcs

import montage_wrapper as montage

def make_template_header(imgdir, c0, pscale=0.2*u.arcsec, xsize=5*u.arcmin,
                         ysize=5*u.arcmin, silent=True, verbose=False):
    '''
    Function to generate ASCII file containing the template header for
    Montage tools

    Parameters
    ----------
    imgdir: string
      Full path to where images are located. Default: None.
      Either sspecify this or fitsfile, but not both

    c0 : `astropy.coordinates` object
      Central coordinate of target

    pscale : astropy.units quantity, optional
      Pixel scale. Default: 0.2 arcsec

    xsize : astropy.units quantity, optional
      Size of desired image in RA. Default: 5 arcmin

    ysize : astropy.units quantity, optional
      Size of desired image in Dec. Default: 5 arcmin

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 2 January 2017
    Modified by Chun Ly, 6 January 2017
     - Bug in sign for CDELT1. Needs to be negative: E to the left.
    '''
    
    if silent == False:
        print '### Begin montage_reproj.make_template_header | '+systime()

    keywords0 = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'CTYPE1',
                 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2', 'CRPIX1',
                 'CRPIX2', 'CROTA2']

    ra0, dec0 = c0.ra.value, c0.dec.value
    nx = np.int(np.ceil((xsize.to(u.arcsec) / pscale.to(u.arcsec)).value))
    ny = np.int(np.ceil((ysize.to(u.arcsec) / pscale.to(u.arcsec)).value))

    values0 = ['T', -64, 2, nx, ny, "'RA---TAN'", "'DEC--TAN'", ra0, dec0,
               -1*pscale.to(u.degree).value, pscale.to(u.degree).value,
               nx/2, ny/2, 0.0]

    txt0 = [a+' = '+str(b) for a,b in zip(keywords0,values0)]
    
    outfile = imgdir + 'template.hdr'
    thefile = open(outfile, 'w')
    for item in txt0: thefile.write(item+'\n')
    thefile.write('END')
    thefile.close()

    if silent == False:
        print '### End montage_reproj.make_template_header | '+systime()

    return outfile
#enddef

def main(c0, fitsfile=None, imgdir=None, out_image=None, xsize=5*u.arcmin,
         ysize=5*u.arcmin, silent=True, verbose=False, catalog='SDSS'):
    '''
    Main function for montage_reproj

    Parameters
    ----------
    fitsfile : string
      Filename of FITS file containing images in different FITS extensions
      Default: None. Either specify this or imgdir, but not both

    imgdir: string
      Full path to where images are located. Default: None.
      Either specify this or fitsfile, but not both

    out_image : string
      Filename of output FITS file. Default is a "final.uncorrected.fits" file
      placed in imgdir or "/var/tmp/montage"

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
    Modified by Chun Ly, 6 January 2017
     - Added out_image keyword option
     - Set ignore=True for call to mProjExec() to avoid crash
     - Delete files to have a fresh start at end of each case
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

    #Mod on 06/01/2017 to handle only multiple images case, output, file deletion
    if n_images > 1:
        t_files = glob.glob(imgdir+'*.fits')
        w0 = WCS(t_files[0])
        pscale0 = np.max(wcs.utils.proj_plane_pixel_scales(w0) * 3600.0) * u.arcsec

        template_header = make_template_header(imgdir, c0, pscale=pscale0,
                                               xsize=xsize, ysize=ysize)
        imgtable = imgdir+'im.tbl'
        montage.mImgtbl(imgdir, imgtable)

        projdir = imgdir+'projdir/'
        if not exists(projdir): os.mkdir(projdir)

        stats_table = projdir + 'stats.tbl'

        # Mod on 06/01/2017 to handle bug with Montage. Ignore stderr
        montage.mProjExec(imgtable, template_header, projdir, stats_table,
                          raw_dir=imgdir, ignore=True)

        proj_imtable = projdir+'im.tbl'
        montage.mImgtbl(projdir, proj_imtable)

        if out_image == None:
            out_image = imgdir+'final.uncorrected.fits'

        montage.mAdd(proj_imtable, template_header, out_image, img_dir=projdir)
        os.system('rm -rf '+imgdir+'im*fits *.tbl *.hdr projdir')

    if silent == False:
        print '### End montage_reproj.main | '+systime()
#enddef
