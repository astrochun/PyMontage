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
                         ysize=5*u.arcmin):
    
    keywords0 = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'CTYPE1',
                 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2', 'CRPIX1',
                 'CRPIX2', 'CROTA2']

    ra0, dec0 = c0.ra.value, c0.dec.value
    nx = np.int(np.ceil((xsize.to(u.arcsec) / pscale.to(u.arcsec)).value))
    ny = np.int(np.ceil((ysize.to(u.arcsec) / pscale.to(u.arcsec)).value))

    values0   = ['T', -64, 2, nx, ny, "'RA---TAN'", "'DEC--TAN'", ra0, dec0,
                 pscale.to(u.degree).value, pscale.to(u.degree).value,
                 nx/2, ny/2, 0.0]

    txt0 = [a+' = '+str(b) for a,b in zip(keywords0,values0)]
    
    outfile = imgdir + 'template.hdr'
    thefile = open(outfile, 'w')
    for item in txt0: thefile.write(item+'\n')
    thefile.write('END')
    thefile.close()

    return outfile
#enddef

def main(c0, fitsfile=None, imgdir=None, xsize=5*u.arcmin, ysize=5*u.arcmin,
         silent=True, verbose=False, catalog='SDSS'):
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

    t_files = glob.glob(imgdir+'*.fits')
    w0 = WCS(t_files[0])
    pscale0 = np.max(wcs.utils.proj_plane_pixel_scales(w0) * 3600.0) * u.arcsec
     
    template_header = make_template_header(imgdir, c0, pscale=pscale0,
                                           xsize=xsize, ysize=ysize)
    imgtable = imgdir+'im.tbl'
    montage.mImgtbl(imgdir, imgtable)

    proj_dir    = imgdir+'projdir/'
    stats_table = proj_dir + 'stats.tbl'

    if not exists(proj_dir): os.mkdir(proj_dir)
    montage.mProjExec(imgtable, template_header, proj_dir, stats_table,
                      raw_dir=imgdir)

    proj_imtable = proj_dir+'im.tbl'
    montage.mImgtbl(projdir, proj_imtable)

    out_image = imgdir+'final.uncorrected.fits'
    print out_image
    montage.mAdd(proj_imtable, template_header, out_image, img_dir=proj_dir)

    if silent == False:
        print '### End montage_reproj.main | '+systime()
#enddef
