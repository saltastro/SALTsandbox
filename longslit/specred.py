
"""
SPECREDUCE

Generdal data reduction script for SALT long slit data.

This includes step that are not yet included in the pipeline 
and can be used for extended reductions of SALT data. 

It does require the pysalt package to be installed 
and up to date.

"""

import os, sys, glob, shutil

import numpy as np
import pyfits
from scipy.ndimage.filters import median_filter


from pyraf import iraf
from iraf import pysalt

from saltobslog import obslog
from saltprepare import saltprepare
from saltbias import saltbias
from saltgain import saltgain
from saltxtalk import saltxtalk
from saltcrclean import saltcrclean
from saltcombine import saltcombine
from saltflat import saltflat
from saltmosaic import saltmosaic
from saltillum import saltillum

from specidentify import specidentify
from specrectify import specrectify
from specsky import skysubtract
from specextract import extract, write_extract
from specsens import specsens
from speccal import speccal

from PySpectrograph.Spectra import findobj

def specred(rawdir, prodir, imreduce=True, specreduce=True, calfile=None, lamp='Ar', automethod='Matchlines', skysection=[800,1000], cleanup=True):
    print rawdir
    print prodir

    #get the name of the files
    infile_list=glob.glob(rawdir+'*.fits')
    infiles=','.join(['%s' % x for x in infile_list])
    

    #get the current date for the files
    obsdate=os.path.basename(infile_list[0])[1:9]
    print obsdate

    #set up some files that will be needed
    logfile='spec'+obsdate+'.log'
    flatimage='FLAT%s.fits' % (obsdate)
    dbfile='spec%s.db' % obsdate

    #create the observation log
    obs_dict=obslog(infile_list)

 
    if imreduce:   
      #prepare the data
      saltprepare(infiles, '', 'p', createvar=False, badpixelimage='', clobber=True, logfile=logfile, verbose=True)

      #bias subtract the data
      saltbias('pP*fits', '', 'b', subover=True, trim=True, subbias=False, masterbias='',  
              median=False, function='polynomial', order=5, rej_lo=3.0, rej_hi=5.0, 
              niter=10, plotover=False, turbo=False, 
              clobber=True, logfile=logfile, verbose=True)

      #gain correct the data
      saltgain('bpP*fits', '', 'g', usedb=False, mult=True, clobber=True, logfile=logfile, verbose=True)

      #cross talk correct the data
      saltxtalk('gbpP*fits', '', 'x', xtalkfile = "", usedb=False, clobber=True, logfile=logfile, verbose=True)

      #cosmic ray clean the data
      #only clean the object data
      for i in range(len(infile_list)):
        if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS'):
          img='xgbp'+os.path.basename(infile_list[i])
          saltcrclean(img, img, '', crtype='edge', thresh=5, mbox=11, bthresh=5.0,
                flux_ratio=0.2, bbox=25, gain=1.0, rdnoise=5.0, fthresh=5.0, bfactor=2,
                gbox=3, maxiter=5, multithread=True,  clobber=True, logfile=logfile, verbose=True)
 
      #flat field correct the data
      flat_imgs=''
      for i in range(len(infile_list)):
        if obs_dict['CCDTYPE'][i].count('FLAT'):
           if flat_imgs: flat_imgs += ','
           flat_imgs += 'xgbp'+os.path.basename(infile_list[i])

      if len(flat_imgs)!=0:
         saltcombine(flat_imgs,flatimage, method='median', reject=None, mask=False,    \
                weight=True, blank=0, scale='average', statsec='[200:300, 600:800]', lthresh=3,    \
                hthresh=3, clobber=True, logfile=logfile, verbose=True)
         saltillum(flatimage, flatimage, '', mbox=11, clobber=True, logfile=logfile, verbose=True)

         saltflat('xgbpP*fits', '', 'f', flatimage, minflat=500, clobber=True, logfile=logfile, verbose=True)
      else:
         flats=None
         imfiles=glob.glob('cxgbpP*fits')
         for f in imfiles:
             shutil.copy(f, 'f'+f)

      #mosaic the data
      geomfile=iraf.osfn("pysalt$data/rss/RSSgeom.dat")
      saltmosaic('fxgbpP*fits', '', 'm', geomfile, interp='linear', cleanup=True, geotran=True, clobber=True, logfile=logfile, verbose=True)

      #clean up the images
      if cleanup:
           for f in glob.glob('p*fits'): os.remove(f)
           for f in glob.glob('bp*fits'): os.remove(f)
           for f in glob.glob('gbp*fits'): os.remove(f)
           for f in glob.glob('xgbp*fits'): os.remove(f)
           for f in glob.glob('fxgbp*fits'): os.remove(f)


    #set up the name of the images
    if specreduce:
       for i in range(len(infile_list)):
           if obs_dict['OBJECT'][i].upper().strip()=='ARC':
               lamp=obs_dict['LAMPID'][i].strip().replace(' ', '')
               arcimage='mfxgbp'+os.path.basename(infile_list[i])
               lampfile=iraf.osfn("pysalt$data/linelists/%s.txt" % lamp)

               specidentify(arcimage, lampfile, dbfile, guesstype='rss', 
                  guessfile='', automethod=automethod,  function='legendre',  order=5, 
                  rstep=100, rstart='middlerow', mdiff=10, thresh=3, niter=5, 
                  inter=True, clobber=True, logfile=logfile, verbose=True)

               specrectify(arcimage, outimages='', outpref='x', solfile=dbfile, caltype='line', 
                   function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
                   blank=0.0, clobber=True, logfile=logfile, verbose=True)
     

    objimages=''
    for i in range(len(infile_list)):
       if obs_dict['CCDTYPE'][i].count('OBJECT') and obs_dict['INSTRUME'][i].count('RSS'):
          if objimages: objimages += ','
          objimages+='mfxgbp'+os.path.basename(infile_list[i])

    if specreduce:
      #run specidentify on the arc files

      specrectify(objimages, outimages='', outpref='x', solfile=dbfile, caltype='line', 
           function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
           blank=0.0, clobber=True, logfile=logfile, verbose=True)


    #create the spectra text files for all of our objects
    spec_list=[]
    for img in objimages.split(','):
       spec_list.extend(createspectra('x'+img, obsdate, smooth=False, skysection=skysection, clobber=True))
    print spec_list
 
    #determine the spectrophotometric standard
    extfile=iraf.osfn('pysalt$data/site/suth_extinct.dat')

    for spec, am, et, pc in spec_list:
        if pc=='CAL_SPST':
           stdstar=spec.split('.')[0]
           print stdstar, am, et
           stdfile=iraf.osfn('pysalt$data/standards/spectroscopic/m%s.dat' % stdstar.lower().replace('-', '_'))
           print stdfile
           ofile=spec.replace('txt', 'sens')
           calfile=ofile #assumes only one observations of a SP standard
           specsens(spec, ofile, stdfile, extfile, airmass=am, exptime=et,
                stdzp=3.68e-20, function='polynomial', order=3, thresh=3, niter=5,
                clobber=True, logfile='salt.log',verbose=True)
    

    for spec, am, et, pc in spec_list:
        if pc!='CAL_SPST':
           ofile=spec.replace('txt', 'spec')
           speccal(spec, ofile, calfile, extfile, airmass=am, exptime=et, 
                  clobber=True, logfile='salt.log',verbose=True)
           #clean up the spectra for bad pixels
           cleanspectra(ofile)

           if spec.count('MCG'): mcg_list.append(ofile)

    #combine the spectra
    #speccombine(mcg_list, obsdate)
    

    #plot the spectra

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
    
    m=(f*0.0)+1
    for i in range(len(m)):
        if f[i]<=0.0:
           x1=int(i-grow)
           x2=int(i+grow)
           m[x1:x2]=0
    m[0]=0
    m[-1]=0

  
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

    
 

def createspectra(img, obsdate, minsize=5, thresh=3, skysection=[800,1000], smooth=False, maskzeros=True, clobber=True):
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

    #find the sections in the images
    section=findobj.findObjects(data, method='median', specaxis=1, minsize=minsize, thresh=thresh, niter=5)
    print section

    #use a region near the center to create they sky
    skysection=findskysection(section, skysection)
    print skysection
 
    #sky subtract the frames
    shdu=skysubtract(hdu, method='normal', section=skysection)
    if os.path.isfile('s'+img): os.remove('s'+img)
    shdu.writeto('s'+img)
 
    spec_list=[]
    #extract the spectra
    #extract the comparison spectrum
    section=findobj.findObjects(shdu[1].data, method='median', specaxis=1, minsize=minsize, thresh=thresh, niter=5)
    print section
    for j in range(len(section)):
        ap_list=extract(shdu, method='normal', section=[section[j]], minsize=minsize, thresh=thresh, convert=True)
        ofile='%s.%s_%i_%i.txt' % (target, obsdate, extract_number(img), j)
        write_extract(ofile, [ap_list[0]], outformat='ascii', clobber=clobber)
        spec_list.append([ofile, airmass, exptime, propcode])

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
   rawdir=sys.argv[1]
   prodir=os.path.curdir+'/'
   agnred(rawdir, prodir)
