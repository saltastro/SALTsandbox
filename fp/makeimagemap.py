import sys
import pyfits
import numpy as np

import pylab as pl

from PySpectrograph.Spectra import detectlines as dl 

def imagemap(img_list, x1, x2, y1, y2, rx,ry):

    nimg=len(img_list)
    dx=x2-x1
    dy=y2-y1
    data=np.zeros((nimg, dy, dx))
    warr=np.zeros(nimg)
    flux=np.zeros(nimg)
 
    rx1=rx-10
    rx2=rx+10
    ry1=ry-10
    ry2=ry+10

    for i,f in enumerate(img_list):
        hdu=pyfits.open(f.strip())
        data[i,:,:]=hdu[1].data[y1:y2,x1:x2]
        warr[i]=hdu[0].header['ET1A']+hdu[0].header['ET1B']*hdu[0].header['ET1Z']
        flux[i]=hdu[1].data[ry1:ry2,rx1:rx2].sum()
        print i, warr[i]
    pl.figure()
    dim=np.zeros((dy,dx))
    for i in range(dx):
     for j in range(dy):
         d=data[:, j,i]/flux
         #dim[j,i]= d.sum() #(warr*d).sum()/ d.sum()
         if d.sum() > 0.025:
            dim[j,i]=dl.centroid(warr, d, kern=[0,-1,0,1,0])
         else:
            dim[j,i]=np.nan
    pl.plot(warr, data[:,300,100])
    pl.plot(warr, data[:,250,150])
    pl.plot(warr, data[:,150,100])
    pl.plot(warr, data[:,100,250])
    d=data[:,300,100]/flux
    print (data[:,100,250]/flux).sum()
    print (warr*(d-d.min())).sum()/(d-d.min()).sum()
    print dl.centroid(warr, d, kern=[0,-1,0,1,0])
    #pl.imshow(dim, vmin=6697, vmax=6722)
    pl.show()

if __name__=='__main__':
   img_list=open(sys.argv[1]).readlines()
   x=int(sys.argv[2])
   y=int(sys.argv[3])
   x1=x-200
   x2=x+200
   y1=y-200
   y2=y+200

   imagemap(img_list, x1, x2, y1, y2, rx=757, ry=263)


    
