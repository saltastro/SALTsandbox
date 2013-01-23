import sys, os
import numpy as np

import pyfits
import numpy

from makeflat import iterstat

def runrandom(img, oimg):

    hdu=pyfits.open(img)
    data=hdu[0].data
    
    #make a new image with filling in the zero pixels with random values
    m,s=iterstat(data[data>0])

    rdata=np.random.normal(loc=m, scale=s,size=data.shape)
    hdu[0].data[data==0]=rdata[data==0]
 
    if os.path.isfile(oimg): os.remove(oimg)
    hdu.writeto(oimg)

if __name__=='__main__':

   infile=open('list.objects').readlines()
   for img in infile:
       img=img.strip()
       img='f'+img
       rimg='r'+img
       runrandom(img, rimg)
       txt=img.replace('fits', 'cat')
       os.system('sex %s -CATALOG_NAME %s' % (rimg, txt))




