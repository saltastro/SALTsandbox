import sys, os
import numpy as np
import pyfits


if __name__=='__main__':


   infile=open('list.objects').readlines()
   for img in infile:
       img=img.strip()
       img='f'+img
       hdu=pyfits.open(img)
       #hdu[0].data[0:250,:]=0  
       hdu[0].data[960:,:]=0

       if os.path.isfile(img): os.remove(img)
       hdu.writeto(img)




