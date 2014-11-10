import sys, os
import numpy as np
from astropy.io import fits 
from astropy import units as u

import ccdproc

def process(ccd, gain, oscan, tsec): 
    """Basic CCD processing for required for data
       
    
    """
    
    #oscan subtract
    ccd = ccdproc.subtract_overscan(ccd, overscan=oscan, median=True)
 
    #gain correct
    ccd = ccdproc.gain_correct(ccd, gain, gain_unit=u.electron/u.adu)

    #trim the image
    ccd = ccdproc.trim_image(ccd, fits_section=tsec)

    return ccd

 

def blueprocess(hdu):
    """Pre-processing necessary for blue files
 
    """
    #first thing first, split the images
    rccd = ccdproc.CCDData(1.0*hdu[0].data[:,0:1050], unit=u.adu)
    lccd = ccdproc.CCDData(1.0*hdu[0].data[:,1050:], unit=u.adu)

    #get the gain values
    gain = hdu[0].header['GAIN'].split()
    #print gain
    
    #process the data
    rccd = process(rccd, float(gain[0]), rccd[:,0:26], '[27:,:]')
    lccd = process(lccd, float(gain[1]), lccd[:,-26:], '[:1024,:]')

     
    #put it back together
    ry,rx = rccd.shape
    ly,lx = lccd.shape
    shape = (ry, rx+lx)
    data = np.zeros(shape)
    data[:,0:rx] = rccd.data
    data[:,rx:rx+lx] = lccd.data

    hdu[0].data = data[::-1,::-1]
    return hdu


if __name__=='__main__':
   infile =sys.argv[1]
   hdu = fits.open(infile)

   hdu = blueprocess(hdu)

   outfile = 'p'+infile
   if os.path.isfile(outfile): os.remove(outfile)
   hdu.writeto(outfile)

     
