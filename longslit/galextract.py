
"""
GALEXTRACT

extract galaxy spectra from image or images and combine them

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter


from pyraf import iraf
from iraf import pysalt

from specsky import skysubtract
from specextract import extract, write_extract
from specsens import specsens
from speccal import speccal

from PySpectrograph.Spectra import findobj

def galextract(imglist, calfile=None):

    
    #set up some files that will be needed
    logfile='specext.log'


    #create the spectra text files for all of our objects
    spec_list=[]
    for img in imglist:
       spec_list.extend(extract_spectra(img, smooth=False, clobber=True))
    print spec_list

    #combine the spectra if desired
 

def speccombine(spec_list, obsdate):
   """Combine N spectra"""

   w1,f1, e1=np.loadtxt(spec_list[0], usecols=(0,1,2), unpack=True)

   w=w1
   f=1.0*f1
   e=e1**2

   for sfile in spec_list[1:]:
      w2,f2, e2=np.loadtxt(sfile, usecols=(0,1,2), unpack=True)
      if2=np.interp(w1, w2, f2)
      ie2=np.interp(w1, w2, e2)
      f2=f2*np.median(f1/if2)
      f+=if2
      e+=ie2**2

   f=f/len(spec_list)
   e=e**0.5/len(spec_list)

   sfile='MCG-6-30-15.%s.spec' % obsdate
   fout=open(sfile, 'w')
   for i in range(len(w)):
           fout.write('%f %e %e\n' % (w[i], f[i], e[i]))
   fout.close()


def cleanspectra(sfile, grow=6):
    """Remove possible bad pixels"""
    try:
        w,f,e=np.loadtxt(sfile, usecols=(0,1,2), unpack=True)
    except:
        w,f=np.loadtxt(sfile, usecols=(0,1), unpack=True)
        e=f*0.0+f.std()

    #set a few unreliable sky lines to zero
    for l in [5577, 6300, 6364]:
        f[abs(w-l)<2]=0
    f[0]=0
    f[-1]=0

    #remove and grow the 
    m=(f*0.0)+1
    for i in range(len(m)):
        if f[i]==0.0:
           x1=int(i-grow)
           x2=int(i+grow)
           m[x1:x2]=0

  
    fout=open(sfile, 'w')
    for i in range(len(w)):
        if m[i]:
           fout.write('%f %e %e\n' % (w[i], f[i], e[i]))
    fout.close()
 
def normalizespectra(sfile, compfile):
    """Normalize spectra by the comparison object"""

    #read in the spectra
    w,f,e=np.loadtxt(sfile, usecols=(0,1,2), unpack=True)
   
    #read in the comparison spectra
    cfile=sfile.replace('MCG-6-30-15', 'COMP')
    print cfile
    wc,fc,ec=np.loadtxt(cfile, usecols=(0,1,2), unpack=True)

    #read in the base star
    ws,fs,es=np.loadtxt(compfile, usecols=(0,1,2), unpack=True)
 
    #calcualte the normalization
    ifc=np.interp(ws, wc, fc) 
    norm=np.median(fs/ifc)
    print norm
    f=norm*f
    e=norm*e

    #write out the result
    fout=open(sfile, 'w')
    for i in range(len(w)):
        fout.write('%f %e %e\n' % (w[i], f[i], e[i]))
    fout.close()

    #copy

    
 

def extract_spectra(img, minsize=5, thresh=3, smooth=False, maskzeros=False, clobber=True):
    """Create a list of spectra for each of the objects in the images"""
    #okay we need to identify the objects for extraction and identify the regions for sky extraction
    #first find the objects in the image
    hdu=pyfits.open(img)
    target=hdu[0].header['OBJECT']
    propcode=hdu[0].header['PROPID']
    airmass=hdu[0].header['AIRMASS']
    exptime=hdu[0].header['EXPTIME']

    if smooth:
       data=smooth_data(hdu[1].data)
    else:
       data=hdu[1].data

    #replace the zeros with the average from the frame
    if maskzeros:
       mean,std=iterstat(data[data>0])
       rdata=np.random.normal(mean, std, size=data.shape)
       print mean, std
       data[data<=0]=rdata[data<=0]

    #use manual intervention to get section
    os.system('ds9 %s &' % img)

    y1=int(raw_input('y1:'))
    y2=int(raw_input('y2:'))
    sy1=int(raw_input('sky y1:'))
    sy2=int(raw_input('sky y2:'))

    ap_list=extract(hdu, method='normal', section=[(y1,y2)], minsize=minsize, thresh=thresh, convert=True)
    sk_list=extract(hdu, method='normal', section=[(sy1,sy2)], minsize=minsize, thresh=thresh, convert=True)
    
    ap_list[0].ldata=ap_list[0].ldata-float(y2-y1)/(sy2-sy1)*sk_list[0].ldata
    
    ofile='%s.%s_%i.txt' % (target, extract_date(img), extract_number(img))
    write_extract(ofile, [ap_list[0]], outformat='ascii', clobber=clobber)

    cleanspectra(ofile)
    
    spec_list=[ofile, airmass, exptime, propcode]

    return spec_list

def smooth_data(data, mbox=25):
    mdata=median_filter(data, size=(mbox, mbox))
    return data-mdata

def find_section(section, y):
    """Find the section closest to y"""
    best_i=-1
    dist=1e5
    for i in range(len(section)):
        d=min(abs(section[i][0]-y), abs(section[i][1]-y))
        if d < dist:
           best_i=i
           dist=d
    return best_i

def extract_number(img):
    """Get the image number only"""
    img=img.split('.fits')
    nimg=int(img[0][-4:])
    return nimg

def extract_date(img):
    """Get the date"""
    img=img.split('.fits')
    obsdate=int(img[0][-8:-4])
    print obsdate
    return obsdate

def iterstat(data, thresh=3, niter=5):
    mean=data.mean()
    std=data.std()
    for i in range(niter):
        mask=(abs(data-mean)<std*thresh)
        mean=data[mask].mean()
        std=data[mask].std()
    return mean, std

    

def findskysection(section, skysection=[800,900], skylimit=100):
    """Look through the section lists and determine a section to measure the sky in

       It should be as close as possible to the center and about 200 pixels wide
    """
    #check to make sure it doesn't overlap any existing spectra
    #and adjust if necessary
    for y1, y2 in section:
        if -30< (skysection[1]-y1)<0:
           skysection[1]=y1-30
        if 0< (skysection[0]-y2)<30:
           skysection[0]=y2+30
    if skysection[1]-skysection[0] < skylimit: print "WARNING SMALL SKY SECTION"
    return skysection
    


if __name__=='__main__':
   galextract(sys.argv[1:])
