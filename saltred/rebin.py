import os,sys 
import pyfits
from pyraf import iraf
from iraf import pysalt

import salttran

def rebin(hdu, xbin, ybin):
    """Given in input SALT image (ie, multi-extension), rebin the data
       and write out the results

       hdu--an input SALT image
       xbin--the factor to rebin the x-data by
       ybin--the factor to rebing the y data by

    """
    for i in range(1,len(hdu)):
          shape=hdu[i].data.shape
          newshape=(ybin*shape[0],xbin*shape[1])
          hdu[i].data=salttran.rebin_factor(hdu[i].data,newshape)

    return hdu
       

if __name__=='__main__':
   inimg=sys.argv[1]
   outimg=sys.argv[2]
   xbin=float(sys.argv[3])
   ybin=float(sys.argv[4])
  
   hdu=pyfits.open(inimg)

   ohdu=rebin(hdu, xbin, ybin)

   if os.path.isfile(outimg):
      os.remove(outimg)

   ohdu.writeto(outimg)

