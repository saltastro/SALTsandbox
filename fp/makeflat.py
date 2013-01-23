import sys, os
import numpy as np
import numpy.ma as ma

import pyfits
from scipy.signal import convolve



def makeflat(data, axc, ayc, arad, order=3, minlevel=10, mbin=20):
    yax,xax=data.shape
    mdata=1.0*data
    mdata=median_radial_mask(mdata, axc, ayc, arad, nbins=50, thresh=3, maskval=0)

    xarr=np.arange(xax)
    flo=minlevel
    fhi=mdata.max()
    for i in range(yax):
       f=mdata[i,:]
       mask=(f>flo)*(f<fhi)
       lmask=(f<flo)
       if mask.sum()>0:
          coef=np.polyfit(xarr[mask], f[mask], order)
          mdata[i,:]=np.polyval(coef, xarr)

    mdata=convolve(mdata, np.ones((mbin,mbin)), mode='same')
    mdata[data<flo]=1
    mdata=data/mdata*mdata[data>flo].mean()
    mdata[data<flo]=0
    return mdata

def median_radial_mask(data, xc=0, yc=0, rmax=500, nbins=50, thresh=3, maskval=0):
    """
    remove pixels which are above the threshhold in the image

    data - 2D image
    xc - x center
    yc - y center
    rmax - maximum radius for the distribution
    nbins - number of bins to use
    maskval - threshold value for including data in the profile
    """

    # calculate the indices from the image
    y, x = np.indices(data.shape)
    radius = np.hypot(x - xc, y - yc)

    #create the radii to sample out to
    hbin=0.5*float(rmax)/nbins
    rpro=np.arange(0,rmax,float(rmax)/nbins)+hbin
    arr=np.zeros(len(rpro))

    #step through each radius and determine the median value for the bin
    for i,r in enumerate(rpro):
        #create the annulus to look at
        mask=(radius>r-hbin)*(radius<r+hbin)
        m,s=iterstat(data[mask], thresh)
        rmask=mask*(abs(data-m)>thresh*s)
        data[rmask]=maskval

    return data

def iterstat(data, thresh=3, niter=5):
    m=np.median(data)
    s=np.std(data)
    for i in range(niter):
        mask=(abs(data-m)<thresh*s)
        m=np.median(data[mask])
        s=np.std(data[mask])
    return m,s


if __name__=='__main__':


   infile=open('list.objects').readlines()
   for img in infile:
       img=img.strip()
       oimg='f'+img
       hdu=pyfits.open(img)
       data=hdu[0].data

       axc=791
       ayc=505
       arad=480

       hdu[0].data=makeflat(data, axc, ayc, arad, order=3, minlevel=10, mbin=20)

       if os.path.isfile(oimg): os.remove(oimg)
       hdu.writeto(oimg)




