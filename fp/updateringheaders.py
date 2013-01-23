import sys, os
import numpy as np
import pyfits

from velpoint import fpfunc

from fit import fit,Parameter

import pylab as pl

"""For each data, look at the information from the calibration file
and update the FITS hearder with the appropriate information
"""

class fitfp:
    def __init__(self, r, r_err,z, t, w, fpcoef, fitcoef=6*[True]):
        self.fpcoef=fpcoef
        self.fitcoef=fitcoef
        for i in range(len(self.fpcoef)):
            if self.fitcoef[i]: self.fpcoef[i]=Parameter(self.fpcoef[i])
        self.r=r
        self.w=w
        self.z=z
        self.t=t-t[0]
        self.r_err=r_err

    def errf(self, w):
        fpcoef=[]
        for i in range(len(self.fpcoef)):
            if self.fitcoef[i]: 
                 fpcoef.append(self.fpcoef[i]())
            else:
                 fpcoef.append(self.fpcoef[i])
     
        return (w-fpfunc(self.z, self.r, self.t, coef=fpcoef))/self.r_err

    def run_fit(self):
        param=[]
        for i in range(len(self.fpcoef)):
            if self.fitcoef[i]: param.append(self.fpcoef[i])
        
        fit(self.errf, param, self.w)

if __name__=='__main__':

   #calculate a list of zeropoints for the data

   #set up the coefficients
   et1A=6741.267818
   et1B=120.6028/1000.0
   fpcoef=[et1A,et1B, 0, 0, 0, 5707.23]
   print fpcoef
   fitcoef=[1, 0, 0, 0, 1, 0]

   #read in the file
   r,r_err,z,t,w=np.loadtxt(sys.argv[1], usecols=(0,1,4,5,6), unpack=True) 
   t=t-t[0]
   wf=fpfunc(z, r, t, coef=fpcoef)
   print wf
   dw=w-wf
   coef=np.polyfit(t, dw, 1)
   print coef

   for i in range(len(t)):
       fpcoef=[et1A+coef[1]+coef[0]*t[i], et1B, 0, 0,0, 5707.23]
       print fpcoef
       wf=fpfunc(z[i], r[i], t[i], coef=fpcoef)
       print w[i], wf

   #add the FP solution to the image headers
   infile=open(sys.argv[2]).readlines()
   for img in infile:
       img=img.strip()
       hdu=pyfits.open(img)
       time=
       fpa=et1A+coef[1]+coef[0]*time
       #header['FPA'] = (, 'Dark Image Subtraction')
