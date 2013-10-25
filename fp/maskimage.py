import sys, os
import numpy as np
import numpy.ma as ma

import pyfits
from scipy.signal import convolve



def maskimage(data, axc, ayc, arad, nbins=50, thresh=3, maskval=0):
    yax,xax=data.shape
    mdata=1.0*data
    mdata=median_radial_mask(mdata, axc, ayc, arad, nbins, thresh, maskval)

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
        if s==0: return np.median(data), np.std(data)
    return m,s


if __name__=='__main__':


   img=sys.argv[1]
   oimg=sys.argv[2]
   axc=float(sys.argv[3])
   ayc=float(sys.argv[4])
   arad=float(sys.argv[5])
   hdu=pyfits.open(img)
   hdu[0].data=maskimage(hdu[0].data, axc, ayc, arad,  nbins=20)

   if os.path.isfile(oimg): os.remove(oimg)
   hdu.writeto(oimg)




