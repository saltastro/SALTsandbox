import sys
import datetime

import pyfits
import numpy as np

import pylab as pl

from matchstars import linearfunc
from makeflat import iterstat

from PySpectrograph.Spectra import detectlines as dl

from fitdata import fitdata

#plot the velocity of a single point from all the images

def fpfunc(z, r, t, coef=None):
    if len(coef)!=6: raise Exception('Not enough FP Coefficients')
    if coef[5]==0:  raise Exception('F must not be zero')
    w=coef[0]+coef[1]*z+coef[2]*z**2+coef[3]*z**3+coef[4]*t
    w=w/(1+(r/coef[5])**2)**0.5
    return w

def qcent(w, f):
    return (w*f).sum()/f.sum()

def invtransform(t_x, t_y, x, y):

    if t_x[1]==0 and t_x[2]==0: 
       return t_x[0]
    elif t_x[2]==0: 
       w=(x-t_x[0])/t_x[1]
    elif t_x[1]==0: 
       w=(y-t_x[0])/t_x[2]
    else:
       w=(x-t_x[0])/t_x[2]-(y-t_y[0])/t_y[2]
       w=w/(t_x[1]/t_x[2]-t_y[1]/t_y[2])
  
    if t_y[1]==0 and t_y[2]==0: 
       return t_y[0]
    elif t_y[2]==0: 
       v=(x-t_y[0])/t_y[1]
    elif t_y[1]==0: 
       v=(y-t_y[0])/t_y[2]
    else:
       v=(x-t_x[0])/t_x[1]-(y-t_y[0])/t_y[1]
       v=v/(t_x[2]/t_x[1]-t_y[2]/t_y[1])
    return w,v

class velpoint:
   def __init__(self, hdulist, axc=0, ayc=0, norm=None, transform=None):
       self.hdulist=hdulist
       self.nframe=len(self.hdulist)


       self.axc=axc
       self.ayc=ayc

       if transform is None:
          self.transform=self.nframe*[[0,1,0,0,0,1]]
          print self.transform
       else:
          self.transform=transform
 
       if norm is None:  
          self.norm=np.ones(self.nframe)
       else:
          self.norm=norm


   def velpoint(self, xr,yr, plot=True):
       wave=np.zeros(self.nframe)
       flux=np.zeros(self.nframe)
       for i in range(self.nframe):
           l=self.transform[i]
           t_x=[float(x) for x in l[0:3]]
           t_y=[float(x) for x in l[3:6]]
 
           #open the data and background subtract it     
           hdu=self.hdulist[i]
           data=hdu[0].data

           #get the header keywords
           et1z=hdu[0].header['ET1Z']
           fpA=hdu[0].header['FPA']
           fpB=hdu[0].header['FPB']
           fpC=hdu[0].header['FPC']
           fpD=hdu[0].header['FPD']
           fpE=hdu[0].header['FPE']
           fpF=hdu[0].header['FPF']
           utctime=hdu[0].header['UTC-OBS']
           utctime=datetime.datetime.strptime(utctime, '%H:%M:%S.%f')  
           fpcoef=[fpA, fpB, fpC, fpD, fpE, fpF] 


           #determine the matching indices
           w,v=invtransform(t_x, t_y, xr, yr)
           radius=((w-self.axc)**2+(v-self.ayc)**2)**0.5
           mid=(int(v), int(w))
           wave[i]=fpfunc(et1z, radius, 0, fpcoef)
           flux[i]=data[mid]*self.norm[i]
 
           #if plot: print xr, yr, w, v, et1z, radius,  wave[i], flux[i], self.norm[i]
       return wave, flux, dl.centroid(wave, flux, kern=[0,-1,0,1,0])
 
 


if __name__=='__main__':
   x=int(sys.argv[1])
   y=int(sys.argv[2])
   norm=np.loadtxt('tmp.normflux', usecols=(1,), unpack=True)
     
   infile=open('tmp.transform').readlines()
   nframe=len(infile)

   transform=[]
   hdulist=[]
   for i, info in enumerate(infile):
       l=info.split()
       img=l[0].strip()
       transform.append(l[1:])
       hdulist.append(pyfits.open(img))   

   vp=velpoint(hdulist, axc=794, ayc=513, norm=norm, transform=transform)
   yi,xi=np.indices((5,5))
   yi=yi+y
   xi=x+xi
   wave=None
   flux=None
   fmax=0
   pl.figure()
   r=[]
   v=[]
   dx=20
   dy=20
   for i in range(x-dx, x+dx):
    for j in range(y-dy, y+dy):
     w,f,c=vp.velpoint(i,j)
     if wave is None: 
         wave=w.copy()
     else:
         wave=np.concatenate((wave,w))
     if flux is None: 
         flux=f.copy()
     else:
         flux=np.concatenate((flux,f))
     #pl.plot(w,f )
     r.append(((i-x)**2+(j-y)**2)**0.5)
     v.append(c)
     if f.max()>fmax: 
        fmax=f.max()
        print i,j,fmax
       
   pl.plot(r,v-np.median(v), ls='', marker='o')   
   pl.show()
