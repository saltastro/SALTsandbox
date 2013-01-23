import sys
import datetime
import time

import pyfits
import numpy as np

import pylab as pl

from matchstars import linearfunc
from makeflat import iterstat

from PySpectrograph.Spectra import detectlines as dl

from fitdata import fitdata

from velpoint import velpoint

if __name__=='__main__':
   xarr=range(600,940)
   yarr=range(315,600)
   #xarr=range(650,730)
   #yarr=range(365,420)
   varr=np.zeros((len(yarr),len(xarr)))
 
   vp=velpoint(axc=794, ayc=513)

   fout=open('alldata.cube', 'w')
   vout=open('vcen.data', 'w')
   tstart=time.time()
   for i in range(len(xarr)):
     print i, (time.time()-tstart)/60.0
     for j in range(len(yarr)):
         wave, flux, vcen=vp.velpoint(xarr[i], yarr[j], plot=False)
         varr[j,i]=vcen
         for k in range(len(wave)):
            fout.write('%5i %5i %6.2f %6.2f\n' % (xarr[i], yarr[j], wave[k], flux[k]))
         vout.write('%5i %5i %6.2f %6.2f %6.2f\n' % (xarr[i], yarr[j], vcen, flux.sum(), flux.std()))
   fout.close()
   vout.close()
