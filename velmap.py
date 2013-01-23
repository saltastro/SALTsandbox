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

from velpoint import velpoint, read_transform

if __name__=='__main__':
   norm=np.loadtxt(sys.argv[1], usecols=(1,), unpack=True)
   hdulist, transform, nframe=read_transform(sys.argv[2])
   axc=int(sys.argv[3])
   ayc=int(sys.argv[4])

   xarr=range(int(sys.argv[5]), int(sys.argv[6]))
   yarr=range(int(sys.argv[7]), int(sys.argv[8]))
   print xarr, yarr
   #xarr=range(650,730)
   #yarr=range(365,420)
   varr=np.zeros((len(yarr),len(xarr)))
 
   vp=velpoint(hdulist, axc=axc, ayc=ayc, norm=norm, transform=transform)
   prefix=sys.argv[1].split('.')[0]
   fout=open(prefix+'_alldata.cube', 'w')
   vout=open(prefix+'_vcen.data', 'w')
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
