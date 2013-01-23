import sys, os
import pyfits
import numpy as np

import scipy.signal as ss

import pylab as pl

from PySpectrograph.Spectra import detectlines as dl 

if __name__=='__main__':
   vcenfile=sys.argv[1]
   x,y,v,f,s=np.loadtxt(vcenfile, unpack=True)
   x1=x.min()
   x2=x.max()
   y1=y.min()
   y2=y.max()
   x=x-x1
   y=y-y1
   xs=x2-x1+1
   ys=y2-y1+1
   print v.size
   print y2-y1, x2-x1

   lcen=np.zeros((ys,xs))
   flux=np.zeros((ys,xs))
   fstd=np.zeros((ys,xs))
 
   for i in range(len(x)):
       lcen[int(y[i])][int(x[i])]=v[i]
       flux[int(y[i])][int(x[i])]=f[i]
       fstd[int(y[i])][int(x[i])]=s[i]
   vcen=2.998e5*(lcen/6564.61-1)

   #vcen=ss.convolve(vcen, np.ones((3,3)))/9
   print flux.mean(), np.median(flux), np.median(fstd)
   mask=(flux<np.median(flux)+5*np.median(fstd))
   mask=(flux**0.5<15)
   #vcen[mask]=np.nan
   #flux[mask]=np.nan
   #vcen[vcen>6724]=np.nan
   print vcen.max(), vcen.min()
   fig=pl.figure()
   ax=fig.add_subplot(111)
   #cim=ax.imshow(vcen)
   #cim=ax.imshow(vcen, vmin=6707, vmax=6722, origin='lower')
   #cim=ax.imshow(np.log10(flux))
   #cim = ax.imshow((flux+500)**0.5)
   #cbar = fig.colorbar(cim) # ticks=[-1, 0, 1])
   #pl.show()


   s=abs(flux)/20.0+20**2 #variance on individual pixel
   s=s**0.5/abs(flux) 
   vcenerr=(s/6564.0*2.998e5)**2+5**2
   vcenerr=vcenerr**0.5
   
   print vcen[140:290,50:210].mean()
   lc=lcen[140:290,50:210]
   sc=s[140:290,50:210]
   print lc.mean(), (lc/sc**2).sum()/(1.0/sc**2).sum(), np.median(sc)

   vcen=vcen.astype(np.float32)
   vcenerr=vcenerr.astype(np.float32)
   hdu = pyfits.PrimaryHDU(vcen[140:290,50:210])
   hdulist = pyfits.HDUList([hdu])
   if os.path.isfile('vcen.fits'): os.remove('vcen.fits')
   hdulist.writeto('vcen.fits')
    
   hdu = pyfits.PrimaryHDU(vcenerr[140:290,50:210])
   hdulist = pyfits.HDUList([hdu])
   if os.path.isfile('vcenerr.fits'): os.remove('vcenerr.fits')
   hdulist.writeto('vcenerr.fits')

   vc=vcen[140:290,50:210] 
   ve=vcenerr[140:290,50:210]
   ym,xm=vc.shape
   print ym, xm
   fout=open('vcen.txt', 'w')
   for j in range((ym)):
     for i in range((xm)):
         fout.write('%i %i %i %i\n' % (i,j,vc[j,i],ve[j,i]))
   fout.close()

   y,x=np.indices(vc.shape)
   r=np.hypot(x-80, y-60)
   f=flux[140:290,50:210]
   ax.plot(r[f>100].ravel(), vc[f>100].ravel(), ls='', marker='o')
   pl.show()
