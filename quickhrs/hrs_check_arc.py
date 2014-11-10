#!/usr/bin/env python
import os
from pylab import *

from PySpectrograph.Spectra import Spectrum

from hrstools import * 
from blueprocess import blueprocess



"""Extract a spectrum from an arc and plot
   up the arc as well
"""
def getfitsspec(dw=0.1, res=0.1):
    wrange = None #[y.min(), y.max()]
    shdu = fits.open('thar.fits')
    ctype1=shdu[0].header['CTYPE1']
    crval1=shdu[0].header['CRVAL1']
    cdelt1=shdu[0].header['CDELT1']
    sf=shdu[0].data
    sw=crval1+cdelt1*np.arange(len(shdu[0].data))
    spec=Spectrum.Spectrum(sw, sf, wrange=wrange, dw=dw, stype='continuum', sigma=res)
    return spec


if __name__=='__main__':

   if len(sys.argv) < 2:
      print """\n
Extract spectrum from arc file and overplot arc file

hrs_check_arc [name] [wavelength]

This will extract a spectrum from an arc file given by
name and plot the order closest to wavelength along
with the actualy spectrum of ThAr around that order

"""

      exit()


   figure()
   hdu = fits.open(sys.argv[1])
   if os.path.basename(sys.argv[1]).startswith('H'):
      hdu = blueprocess(hdu)
   wavelength = float(sys.argv[2])
 
   obsmode = hdu[0].header['OBSMODE']
   sp, order_file, wave_file, dy = get_info(wavelength, obsmode)

   orders = fits.open(order_file)[0].data
   arc = 1.0*hdu[0].data
   o=find_order(sp, wavelength)
   print "Order: ", o
   sp.grating.order = o
   data = hdu[0].data
   #todo replace with arc image
   w,flux = get_order(data, orders, arc, sp, o, dy=dy)
   plot(w, flux)
   spec = getfitsspec()
   mask = (spec.wavelength > w.min()) * (spec.wavelength < w.max())
   plot(spec.wavelength[mask], spec.flux[mask]*flux.max()/spec.flux[mask].max(), color='black')
   show()


