import sys, os
import numpy as np
import pyfits
from makeflat import iterstat


if __name__=='__main__':


   infile=open(sys.argv[1]).readlines()
   for img in infile:
       img=img.strip()
       hdu=pyfits.open(img)

       data=hdu[0].data
       m,s=iterstat(data[data>0])
       hdu[0].data[data>0]=data[data>0]-m
       
       oimg=img
       if os.path.isfile(oimg): os.remove(oimg)
       hdu.writeto(oimg)




