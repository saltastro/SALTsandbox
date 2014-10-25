#! /usr/bin/env python

#########################################################################################################
#  REDUCE SALT SPECTRA with IRAF tasks:                                                                 #
#  Info from a wiki Alexei wrote: https://sciencewiki.salt.ac.za/index.php/Long_Slit_Reduction_Recipe   #
#########################################################################################################

# IMPORT all IRAF tasks required for this script
import sys
import os
import numpy
import pyraf
from pyraf import iraf
from pyraf.iraf import images
from pyraf.iraf import immatch
from pyraf.iraf import imutil
from pyraf.iraf import pysalt
from pyraf.iraf import saltred

os.system('mkdir SlitViews')
os.system('mv mbxgpS*fits SlitViews')
os.system('rm Reduced/*fits')
os.system('rm Reduced/*spectra*.txt')
os.system('ls mbxgpP*fits > list')
os.system('mkdir Reduced')
os.system('dfits mbxgpP*fits | fitsort object camang grating ccdsum date-obs time-obs > summary')
filenames = numpy.loadtxt('list',dtype='str')

#  Create MASTER script to handle CALIBRATIONS

calibrate = open('calibrate.py','w')
print >> calibrate, '#! /usr/bin/env python'
print >> calibrate, 'import pyraf'
print >> calibrate,  'from pyraf import iraf'
print >> calibrate,  'from pyraf.iraf import noao'
print >> calibrate,  'from pyraf.iraf import onedspec'
print >> calibrate,  'from pyraf.iraf import twodspec'
print >> calibrate,  'from pyraf.iraf import longslit'

#  Create MASTER script to handle REDUCTIONS

reduction = open('reduction.py','w')
print >> reduction, '#! /usr/bin/env python'
print >> reduction, 'import sys'
print >> reduction, 'import os'
print >> reduction, 'import pyraf'
print >> reduction,  'from pyraf import iraf'
print >> reduction,  'from pyraf.iraf import noao'
print >> reduction,  'from pyraf.iraf import onedspec'
print >> reduction,  'from pyraf.iraf import twodspec'
print >> reduction,  'from pyraf.iraf import longslit'
print >> reduction,  'from pyraf.iraf import images'
print >> reduction,  'from pyraf.iraf import imutil'
print >> reduction,  'from pyraf.iraf import pysalt'
print >> reduction,  'from pyraf.iraf import saltred'
#print >> reduction,  "os.system('rm vbctf*fits bctf*fits ctf*fits tf*fits f*fits')"
print >> reduction,  "os.system('rm *.spectra*.fits fc.Arc*.*')"
print >> reduction,  "os.system('rm *spectra*.txt')"

# Create lists of TARGET/OBJECT, ARC and FLAT
targetlist = open('targetlist','w')
arclist = open('arclist','w')
flatlist = open('flatlist','w')

# Set start values
arccount = 0
targetcount = 0
flatcount = 0
targetfiles = []
filenumber = []
prevcamangle = ''

# Determine how to handle each RSS files in the directory
for i in range(len(filenames)):
   obstype = imutil.hselect(images = filenames[i]+'[0]', fields = 'CCDTYPE', expr = 'yes',Stdout=1)
   obsmode = imutil.hselect(images = filenames[i]+'[0]', fields = 'OBSMODE', expr = 'yes',Stdout=1)
   grating = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'GRATING', expr = 'yes',Stdout=1))
   grangle = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'GR-ANGLE', expr = 'yes',Stdout=1))
   arangle = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'AR-ANGLE', expr = 'yes',Stdout=1))
   camangle = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'CAMANG', expr = 'yes',Stdout=1))
   ccdsum = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'CCDSUM', expr = 'yes',Stdout=1))
   utctime = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'TIME-OBS', expr = 'yes',Stdout=1))
   utcdate = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'DATE-OBS', expr = 'yes',Stdout=1))
   print filenames[i]
   gratingvalues = 'GRATING = '+grating+ '\t'+'GR-ANGLE = '+grangle+'\t'+'AR-ANGLE = '+arangle+'\t'+'CAMANG = '+camangle+'\t'+'CCDSUM = '+ccdsum
   print gratingvalues
   # ARC file:			CCDTYPE = 'ARC               ' 
   if obstype.count('ARC') > 0:
      # Strip multi header:
      imutil.imcopy(input=filenames[i]+'[SCI,inherit]', output='Reduced/Arc'+str(arccount+1)+'.fits',verbose='yes')
      print filenames[i]+' is an ARC'
      arctype = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'LAMPID', expr = 'yes',Stdout=1))
      arctype = arctype.strip('[').strip(']').replace('"','').replace("'",'').replace(' ','')
      print 'LAMPID  = '+arctype
      if arctype == 'Ne':
         arctype = 'NeAr'
      if arctype == 'Ar':
         if float(grating.replace('PG','')) >= 1800:
            arctype = 'Argon_hires'
         else:
            arctype = 'Argon_lores' 
      # Identify known spectral lines in the ARC
     #Alexei wiki: identify images=refspec.fits coordlist=cuar.dat function=chebyshev order=5 fwidth=6. cradius=6
      print >> calibrate, "onedspec.identify(images='cr.Arc"+str(arccount+1)+".fits', coordlist='"+arctype+".txt', function='chebyshev', order=5, fwidth=6, cradius=6, autowrite='yes')"
      # Re-Identify known spectral lines in the ARC
     #Alexei wiki: reidentify reference=refspec.fits images=refspec.fits interactive=no newaps=yes override=no refit=yes nlost=20 coordlist=cuar.dat verbose=yes
      print >> calibrate, "onedspec.reidentify(reference='cr.Arc"+str(arccount+1)+".fits', images='cr.Arc"+str(arccount+1)+".fits', interactive='no', newaps='yes', override='no', refit='yes', nlost=20, coordlist='"+arctype+".txt', verbose='yes')"
      # Obtain theoretical ARC
      os.system('wget http://pysalt.salt.ac.za/lineatlas/'+arctype+'.txt')
      # Fit coordinates of the ARC to wavelengths
     #Alexei wiki: fitcoords images=refspec interactive=yes combine=no functio=legendre xorder=5 yorder=3
      print >> calibrate, "longslit.fitcoords(images='cr.Arc"+str(arccount+1)+"', interac='yes', combine='no', function='legendre', xorder=5, yorder=3)"
      print >> arclist, 'Reduced/Arc'+str(arccount+1)+'.fits'
      # Transform the Arc itself to test the accuracy of curvature correction
      print >> reduction, "longslit.transform(input='cr.Arc"+str(arccount+1)+".fits', output='fc.Arc"+str(arccount+1)+".fits', fitnames='cr.Arc"+str(arccount+1)+"', interptype='linear')"
      arccount = arccount + 1
   # TARGET/OBJECT file:	CCDTYPE = 'OBJECT               '       
   if obstype.count('OBJECT') > 0 and obsmode.count('SPECTROSCOPY') > 0 :
      filenumber.append(filenames[i].strip('.fits').strip('mbxgpP'))
      if targetcount == 0:
         os.system("ds9 "+filenames[i]+"&")
         targetbox = raw_input("Enter y-range for extraction of target in the format y1:y2 : ")
         targetmin = float(targetbox.split(':')[0])
         targetmax = float(targetbox.split(':')[1])
         boxmin = int(targetmin-40)
         boxmax = int(targetmax+40)
         trimbox = str(boxmin)+':'+str(boxmax)
      else:
         if float(filenumber[-1]) - float(filenumber[-2]) > 1:
            print '############################################################################'
            print "  OBJECT files don't follow chronologically. You must reduce it seperately!"
            print '############################################################################'
      if camangle != prevcamangle and prevcamangle != '':
         print '######################################################################'
         print "  Different setting was also observed. You must reduce it seperately!"
         print '######################################################################'
      prevcamangle = camangle
      # Strip multi header:
      imutil.imcopy(input=filenames[i]+'[SCI,inherit]', output='Reduced/spectra'+str(targetcount+1)+'.fits',verbose='yes')
      # TRANSFORM science files
     #Alexei wiki: transform input=stdstar.fits output=stdstarw.fits fitnames=refspec interptype=linear flux=yes blank=INDEF x1=INDEF x2=INDEF dx=INDEF y1=INDEF y2=INDEF dy=INDEF
      print >> reduction, "longslit.transform(input='spectra"+str(targetcount+1)+".fits', output='f.spectra"+str(targetcount+1)+".fits', fitnames='cr.Arc1', interptype='linear')"
      # TRIM science files to appropriate portion only (relevant line and surrounds)
      print >> reduction, "imutil.imcopy(input='f.spectra"+str(targetcount+1)+".fits[*,"+trimbox+"]', output='tf.spectra"+str(targetcount+1)+".fits', verbose='yes')"
      print >> targetlist, 'Reduced/spectra'+str(targetcount+1)+'.fits'
      #print >> targetlist, '# Original filename:', filenames[i], utcdate, utctime
      #print >> targetlist, '#', gratingvalues
      targetcount = targetcount + 1
      targetfiles.append(filenames[i])
   # FLAT file:			CCDTYPE = 'FLAT              '
   if obstype.count('FLAT') > 0:
      # Strip multi header:
      flatfield = 'Reduced/flat'+str(flatcount+1)+'.fits'
      imutil.imcopy(input=filenames[i]+'[SCI,inherit]', output=flatfield,verbose='yes')
      averageflux = imutil.imstat(images=flatfield, mode="h", fields='image,mean,stddev',Stdout=1)[1].split()[-2]
      imutil.imarith(operand1=flatfield, op="/", operand2=str(averageflux), result=flatfield, title="Normalized FLAT", mode="h")
      print filenames[i]+' is a FLAT'
      print >> flatlist, flatfield
      #print >> flatlist, '# Original filename:', filenames[i], gratingvalues
      flatcount = flatcount + 1

if flatcount > 0:
   # Make MASTER flat
   flatlist.close()
   immatch.imcombine.input = '@flatlist'
   immatch.imcombine.output = 'MasterFlat.fits'                               
   immatch.imcombine.combine = "median"
   immatch.imcombine.reject = "none"
   immatch.imcombine.mode = "h"
   immatch.imcombine()

   # Apply illumination correction to master flat
   saltred.saltillum(images='MasterFlat.fits', outpref="", outimages="iMasterFlat.fits", mbox=11, clobber='yes')

   # Flat-field the target frames
   targetlist.close()
   imutil.imarith.operand1 = '@targetlist'
   imutil.imarith.op = "/"
   imutil.imarith.operand2 = 'iMasterFlat.fits'
   imutil.imarith.result = '@targetlist'
   imutil.imarith.mode = "h"
   imutil.imarith()

# Remove Cosmic Rays from Arcs:
saltred.saltcrclean(images='Reduced/*Arc*.fits',outpref='Reduced/cr.',crtype='fast',flux_ratio=0.25)

os.system('rm *.txt.*')
os.system('mv *.txt Reduced/')
os.system('chmod a+x calibrate.py')
os.system('mv calibrate.py Reduced/')
os.system('chmod a+x reduction.py')
os.system('mv reduction.py Reduced/')

# Cosmic ray cleaning : FAST 
#print >> reduction, 'fluxratio = raw_input("Enter flux ratio for Cosmic Ray cleaning task (default: 0.25) : ")'
#print >> reduction, "if fluxratio == '':"
#print >> reduction, '   fluxratio = 0.25'
#print >> reduction, 'else:'
#print >> reduction, '   fluxratio = float(fluxratio)'
#print >> reduction, "saltred.saltcrclean(images='tf.*.fits', outpref='c', crtype='fast', flux_ratio=str(fluxratio))"
#print >> reduction, "saltred.saltcrclean(images='tf.*.fits', outimages='', outpref='c', crtype='fast', thresh=5.0, mbox=5, bthresh=3.0, flux_ratio=str(fluxratio), bbox=11, gain=1.0, rdnoise=5.0, fthresh=5.0, bfactor=2, gbox=3, maxiter=5, multithread='no', clobber='no', logfile='salt.log', verbose='yes')"
# Cosmic ray cleaning : EDGE
#  Parameters suggested in: https://sciencewiki.salt.ac.za/index.php/SALT_Longslit_Reduction_Tutorial
print >> reduction, "saltred.saltcrclean(images='tf.*.fits', outpref='c', crtype='edge', bthresh=5.0, thresh=10.0, maxiter=5, mbox=25, flux_ratio=0.2, bbox=25, gain=1.0, rdnoise=5.0, fthresh=5.0, bfactor=2, gbox=3, multithread='no', clobber='yes', logfile='salt.log', verbose='yes')"

# Create residual images for inspection of Cosmic Rays removed
for j in range(len(targetfiles)):
   print >> reduction, "imutil.imarith(operand1='tf.spectra"+str(j+1)+".fits', op='-', operand2='ctf.spectra"+str(j+1)+".fits', result='cr.residual"+str(j+1)+".fits', title='residual', mode='h')"

# BACKGROUND subtraction
for j in range(len(targetfiles)):
   #Alexei wiki: background input=stdstarw.fits output=stdstarws.fits axis=2 interactive=no naverage=1 function=chebyshev order=5 low_rej=2 high_rej=1.5 niterate=5 grow=0
   print >> reduction, "longslit.background(input='ctf.spectra"+str(j+1)+".fits', output='bctf.spectra"+str(j+1)+".fits', axis=2, interactive='no', naverage=1, function='chebyshev', order=5, low_rej=2, high_rej=1.5, niterate=5, grow=0)"

# CALIBRATE using a spectrophotometric standard (in data folders: CAL_SPST)
#for c in range(len(targetfiles)):
#   print >> reduction, "onedspec.calibrate(input='bctf.spectra"+str(c+1)+".fits', output='bctf0.spectra"+str(c+1)+".fits', sensitivity=sensitivity.fits)"

# Set APERTURES for extraction (looks for largest peak and will fail to determine correctly if peak is saturated)
print >> reduction, 'from pyraf.iraf import apextract'
#Alexei wiki: apall input=stdstarws.fits format=onedspec interactive=no nfind=1 llimit=-40 ulimit=40 t_order=10 ylevel=INDEF
print >> reduction, "apextract.apall(input='bctf.*.fits', interactive='no', nfind=1, upper=10, lower=-10, b_sample = '-30:-20,20:30')"
print >> reduction, "apextract.apall(input='ctf.*.fits', interactive='no', nfind=1, upper=10, lower=-10, b_sample = '-30:-20,20:30')"
#interactive='yes', find='yes', recenter='yes', resize='yes', edit='yes', trace='yes', fittrace='yes', extract='yes', extras='yes', review='yes'
# Extract the 2D observed Arc
print >> reduction, "apextract.apall(input='fc.Arc*.fits', interactive='no', nfind=1, upper=10, lower=-10, b_sample = '-50:-40,40:50')"

# Do Vacuum/Air correction 
for v in range(len(targetfiles)):
   # Background-subtracted images
   print >> reduction, "onedspec.disptrans(input='bctf.spectra"+str(v+1)+".ms.fits', output='vbctf.spectra"+str(v+1)+".ms.fits', units='angs', air='air2vac')"
   # Without background-subtraction
   print >> reduction, "onedspec.disptrans(input='ctf.spectra"+str(v+1)+".ms.fits', output='vctf.spectra"+str(v+1)+".ms.fits', units='angs', air='air2vac')"
for vv in range(arccount):
   # Arc 
   print >> reduction, "onedspec.disptrans(input='fc.Arc"+str(vv+1)+".ms.fits', output='vfc.Arc"+str(vv+1)+".ms.fits', units='angs', air='air2vac')"

# Create 1D spectra in ASCII format from 2D spectra
for k in range(len(targetfiles)):
   # Use Vacuum-corrected 2D spectra
   print >> reduction, "onedspec.wspectext(input='vbctf.spectra"+str(k+1)+".ms.fits', output='vbctf.spectra"+str(k+1)+".txt', header='no', wformat="+'"'+"%0.2f"+'"'+")"
   # Use non Vacuum-corrected 2D spectra
   print >> reduction, "onedspec.wspectext(input='bctf.spectra"+str(k+1)+".ms.fits', output='bctf.spectra"+str(k+1)+".txt', header='no', wformat="+'"'+"%0.2f"+'"'+")"
   # Use Vacuum-corrected 2D spectra Without background-subtraction (to assess skylines to check calibration accuracy)
   print >> reduction, "onedspec.wspectext(input='vctf.spectra"+str(k+1)+".ms.fits', output='vctf.spectra"+str(k+1)+".txt', header='no', wformat="+'"'+"%0.2f"+'"'+")"
   # Use 2D spectra Without vacuum-correction or background-subtraction (to check the necessity of the air2vac correction using skylines)
   print >> reduction, "onedspec.wspectext(input='ctf.spectra"+str(k+1)+".ms.fits', output='ctf.spectra"+str(k+1)+".txt', header='no', wformat="+'"'+"%0.2f"+'"'+")"
for kk in range(arccount):
   # Extract the 1D observed Arc (without vacuum correction)
   print >> reduction, "onedspec.wspectext(input='fc.Arc"+str(kk+1)+".ms.fits', output='fc.Arc"+str(kk+1)+".txt', header='no', wformat="+'"'+"%0.2f"+'"'+")"
   # Extract the 1D observed Arc (with vacuum correction)
   print >> reduction, "onedspec.wspectext(input='vfc.Arc"+str(kk+1)+".ms.fits', output='vfc.Arc"+str(kk+1)+".txt', header='no', wformat="+'"'+"%0.2f"+'"'+")"

###################################################
#  PLOT the 1D spectra together (trailed spectra) #
###################################################

# Prepare GNUPLOT script for the trailed spectra
plots = open('Reduced/plots','w')

dateobs = str(imutil.hselect(images = filenames[0]+'[0]', fields = 'DATE-OBS', expr = 'yes',Stdout=1)).strip('[').strip(']').replace("'",'')

print >> plots, "set output 'Spectra.eps'"
print >> plots, "set terminal postscript portrait"+'"Helvetica" 10'
print >> plots, "set multiplot layout "+str(targetcount+1)+",1 title 'RSS long-slit spectra on "+dateobs+"'"
print >> plots, "set tmargin 0.25"
print >> plots, "set bmargin 0"
print >> plots, "set lmargin at screen 0.15"
print >> plots, "set rmargin at screen 0.98"
print >> plots, "set pointsize 0.3"
print >> plots, "unset xtics"
print >> plots, "set ylabel '  '"
print >> plots, "set label 'Flux' at screen 0.02,0.5 rotate by 90"
print >> plots, "set ytics 0,200,10000"
timeobs = str(imutil.hselect(images = targetfiles[0]+'[0]', fields = 'TIME-OBS', expr = 'yes',Stdout=1)).strip('[').strip(']').replace("'",'')
timeobs = 'UTC '+timeobs
print >> plots, "plot [] [0:] 'vbctf.spectra"+str(1)+".txt' with lines title '"+timeobs+"    ("+targetfiles[0]+")'"
print >> plots, "unset label"
print >> plots, "unset title"
print >> plots, "set tmargin 0"
for l in range(1,len(targetfiles)-1):
   timeobs = str(imutil.hselect(images = targetfiles[l]+'[0]', fields = 'TIME-OBS', expr = 'yes',Stdout=1)).strip('[').strip(']').replace("'",'')
   print >> plots, "plot [] [0:] 'vbctf.spectra"+str(l+1)+".txt' with lines title '"+timeobs+"    ("+targetfiles[l]+")'"
print >> plots, "set xlabel 'Wavelength [Angstrom]'"
print >> plots, "set xtics 0,200,10000"
timeobs = str(imutil.hselect(images = targetfiles[l+1]+'[0]', fields = 'TIME-OBS', expr = 'yes',Stdout=1)).strip('[').strip(']').replace("'",'')
print >> plots, "plot [] [0:] 'vbctf.spectra"+str(len(targetfiles))+".txt' with lines title '"+timeobs+"    ("+targetfiles[l+1]+")'"
print >> plots, "unset multiplot"
print >> plots, "exit"
plots.close()

print >> reduction, "os.system('gnuplot<plots')"
