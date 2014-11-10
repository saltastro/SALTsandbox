#!/usr/bin/env python
import os, sys
from pylab import *

from PySpectrograph.Spectra import Spectrum

from hrstools import * 
from blueprocess import blueprocess


"""Extract a spectrum from an file
"""


if __name__=='__main__':

   if len(sys.argv) < 2:
      print """\n
Extract spectrum from arc file and overplot arc file

hrs_extract_spectrum [name] [wavelength]

This will extract a spectrum from file given by
name and plot the order closest to wavelength.  The
infile should be the raw file for the arm.  Any
processing will be done by the script

"""

      exit()


   figure()
   hdu = fits.open(sys.argv[1])
   if os.path.basename(sys.argv[1]).startswith('H'):
      hdu = blueprocess(hdu)
   data = hdu[0].data

         
   wavelength = float(sys.argv[2])

 
   obsmode = hdu[0].header['OBSMODE']
   sp, order_file, wave_file, dy = get_info(wavelength, obsmode)

   orders = fits.open(order_file)[0].data
   arc = 1.0*hdu[0].data
   o=find_order(sp, wavelength)
   print "Order: ", o
   sp.grating.order = o

   w,flux = get_order(data, orders, arc, sp, o, dy=dy)
   plot(w, flux)
   show()


