import sys

import pylab as pl
import numpy as np

import PySpectrograph
from PySpectrograph.Models import RSSModel
from PySpectrograph.Spectra import Spectrum

def plotarcspectra(arclist='Ne.txt', grating='PG0900', gratang=15.0, camang=30.0, slit=1.5, xbin=2, ybin=2):
   """Plot an Arc Line list for a given RSS set up 
      arclist--an arc line list in the format of 'wavelength flux' for arc lines
      grating--name of the grating
      gratang--angle of the grating
      camang--angle of the camera (articulation angle)
      slit--slit width in arcseconds
      xbin--xbinning
      ybin--ybinning

   """
   rss=RSSModel.RSSModel(grating_name=grating, gratang=gratang, camang=camang, 
                      slit=slit, xbin=xbin, ybin=ybin)


   #print out some basic statistics
   print 1e7*rss.calc_bluewavelength(), 1e7*rss.calc_centralwavelength(), 1e7*rss.calc_redwavelength()
   R=rss.calc_resolution(rss.calc_centralwavelength(), rss.alpha(), -rss.beta())
   res=1e7*rss.calc_resolelement(rss.alpha(), -rss.beta())
   print R, res

   #set up the detector
   ycen=rss.detector.get_ypixcenter()
   d_arr=rss.detector.make_detector()[ycen,:] #creates 1-D detector map
   xarr=np.arange(len(d_arr))
   w=1e7*rss.get_wavelength(xarr)

   #set up the artificial spectrum
   sw,sf=pl.loadtxt(arclist, usecols=(0,1), unpack=True)
   wrange=[1e7*rss.calc_bluewavelength(), 1e7*rss.calc_redwavelength()]
   spec=Spectrum.Spectrum(sw, sf, wrange=wrange, dw=res/10, stype='line', sigma=res)

   #interpolate it over the same range as the detector
   spec.interp(w)
   


   #plot it
   pl.figure()
   pl.plot(spec.wavelength, d_arr*((spec.flux)/spec.flux.max()))
   pl.plot(spec.wavelength, d_arr*0.1)

   print pl.gca().get_ylim()
   for i,f in zip(sw,sf):
       if i>spec.wavelength.min() and i<spec.wavelength.max():
          print i, f, sf.max(), spec.flux.max()
          #pl.text(i, f/sf.max(), i, fontsize='large', rotation=90)
          y=max(0.1, f/sf.max())
          pl.text(i, y, i,  rotation=90)
   
   pl.show()


if __name__=='__main__':
   """Example of running the task
   python plotarcspectra.py [line list] [grating name] [grating angle] [camera angle] [slit] [xbin] [ybin]

   Here is an example of running it:
   python plotarcspectra.py Ne.txt PG0900  15.0 30.0 1.5 2 2

   """
   #create the spectrograph model
   arclist=sys.argv[1]
   grating=sys.argv[2]
   gratang=float(sys.argv[3])
   camang=float(sys.argv[4])
   slit=float(sys.argv[5])
   xbin=int(sys.argv[6])
   ybin=int(sys.argv[7])
   plotarcspectra(arclist, grating, gratang, camang, slit, xbin, ybin)
