#! /usr/bin/env python

##################################################################################################
#  REDUCE SALT SPECTRA with PySALT:                                                              #
#  As suggested in: https://sciencewiki.salt.ac.za/index.php/SALT_Longslit_Reduction_Tutorial    #
##################################################################################################

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
print >> calibrate,  'from pyraf.iraf import pysalt'
print >> calibrate,  'from pyraf.iraf import saltred'
print >> calibrate,  'from pyraf.iraf import saltspec'

#  Create MASTER script to handle REDUCTIONS

reduction = open('reduction.py','w')
print >> reduction, '#! /usr/bin/env python'
print >> reduction, 'import sys'
print >> reduction, 'import os'
print >> reduction, 'import pyraf'
print >> reduction,  'from pyraf import iraf'
print >> reduction,  'from pyraf.iraf import images'
print >> reduction,  'from pyraf.iraf import imutil'
print >> reduction,  'from pyraf.iraf import pysalt'
print >> reduction,  'from pyraf.iraf import saltred'
print >> reduction,  'from pyraf.iraf import saltspec'
print >> reduction,  "os.system('rm bctf*fits ctf*fits tf*fits fc*fits')"
print >> reduction,  "os.system('rm *spectra*.txt')"

# Create lists of TARGET/OBJECT, ARC and FLAT
targetlist = open('targetlist','w')
flatfieldedlist = open('flatfieldedlist','w')
arclist = open('arclist','w')
flatlist = open('flatlist','w')
inputs = open('inputs','w')

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
      # Rename files:
      os.system('cp '+filenames[i]+' Reduced/Arc'+str(arccount+1)+'.fits')
      print filenames[i]+' is an ARC'
      arctype = str(imutil.hselect(images = filenames[i]+'[0]', fields = 'LAMPID', expr = 'yes',Stdout=1))
      arctype = arctype.strip('[').strip(']').replace('"','').replace("'",'').replace(' ','')
      grating = grating.strip('[').strip(']').replace('"','').replace("'",'').replace(' ','')
      print 'LAMPID  = '+arctype
      if arctype == 'Ne':
         arctype = 'NeAr'
      if arctype == 'Ar':
         if float(grating.replace('PG','')) >= 1800:
            arctype = 'Argon_hires'
         else:
            arctype = 'Argon_lores'   
      # Identify known spectral lines in the ARC and Fit coordinates of the ARC to wavelengths
      print >> calibrate, "saltspec.specidentify(images='cr.Arc"+str(arccount+1)+".fits', linelist='"+arctype+".txt', outfile='db"+str(arccount+1)+".sol', guesstype='rss', guessfile='', automethod='Matchlines', function='legendre', order=3, rstep=200, rstart='middlerow', mdiff=20, thresh=3, niter=5, startext=0, inter='yes', clobber='no', logfile='salt.log', verbose='yes')"     
      # Obtain theoretical ARC
      os.system('wget http://pysalt.salt.ac.za/lineatlas/'+arctype+'.txt')
      print >> arclist, 'Reduced/Arc'+str(arccount+1)+'.fits'
      print >> reduction, "saltspec.specrectify(images='cr.Arc"+str(arccount+1)+".fits', outimages='', outpref='t', solfile='db"+str(arccount+1)+".sol', caltype='line', function='legendre', order=3, inttype='interp', w1='None', w2='None', dw='None', nw='None', blank=0.0, clobber='yes', logfile='salt.log', verbose='yes')"
      arccount = arccount + 1
   # TARGET/OBJECT file:	CCDTYPE = 'OBJECT               '       
   if obstype.count('OBJECT') > 0 and obsmode.count('SPECTROSCOPY') > 0 :
      filenumber.append(filenames[i].strip('.fits').strip('mbxgpP'))
      if targetcount == 0:
         os.system("ds9 "+filenames[i]+"&")
         targetbox = raw_input("Enter y-range for extraction of TARGET in the format y1:y2 : ")
         targetmin = float(targetbox.split(':')[0])
         targetmax = float(targetbox.split(':')[1])
         skybox =  raw_input("Enter y-range for extraction of SKY in the format y1:y2 : ")
         skymin = float(skybox.split(':')[0])
         skymax = float(skybox.split(':')[1])
         boxmin = min(targetmin,skymin)-40
         boxmax = max(targetmax,skymax)+40
         trimbox = str(boxmin)+':'+str(boxmax)
         print >> inputs, targetbox
         print >> inputs, skybox
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
      # Rename files:
      os.system('cp '+filenames[i]+' Reduced/spectra'+str(targetcount+1)+'.fits')
      # Strip multi header:
      #imutil.imcopy(input=filenames[i]+'[SCI,inherit]', output='Reduced/spectra'+str(targetcount+1)+'.fits',verbose='yes')
      # TRANSFORM science files
      print >> reduction, "saltspec.specrectify(images='f.spectra"+str(targetcount+1)+".fits', outimages='', outpref='t', solfile='db1.sol', caltype='line', function='legendre', order=3, inttype='interp', w1='None', w2='None', dw='None', nw='None', blank=0.0, clobber='yes', logfile='salt.log', verbose='yes')"
      print >> targetlist, 'Reduced/spectra'+str(targetcount+1)+'.fits'
      print >> flatfieldedlist, 'Reduced/f.spectra'+str(targetcount+1)+'.fits'
      #print >> targetlist, '# Original filename:', filenames[i], utcdate, utctime
      #print >> targetlist, '#', gratingvalues
      targetcount = targetcount + 1
      targetfiles.append(filenames[i])
   # FLAT file:			CCDTYPE = 'FLAT              '
   if obstype.count('FLAT') > 0:
      # Rename files:
      flatfield = 'Reduced/flat'+str(flatcount+1)+'.fits'
      os.system('cp '+filenames[i]+' '+flatfield)
      print filenames[i]+' is a FLAT'
      print >> flatlist, flatfield
      #print >> flatlist, '# Original filename:', filenames[i], gratingvalues
      flatcount = flatcount + 1

inputs.close()
os.system('cp inputs oldinputs')

if flatcount > 0:
   # Make MASTER flat (saltcombine automatically normalizes the flats)
   flatlist.close()
   saltred.saltcombine.images = '@flatlist'
   saltred.saltcombine.outimage = 'MasterFlat.fits'                               
   saltred.saltcombine.combine = "median"
   saltred.saltcombine.reject = "None"
   saltred.saltcombine.mode = "h"
   saltred.saltcombine()

   # Apply illumination correction to master flat
   saltred.saltillum(images='MasterFlat.fits', outpref="", outimages="iMasterFlat.fits", mbox=11, clobber='yes')

   # Flat-field the target frames
   targetlist.close()
   flatfieldedlist.close()
   saltred.saltarith.operand1 = '@targetlist'
   saltred.saltarith.op = "/"
   saltred.saltarith.operand2 = 'iMasterFlat.fits'
   saltred.saltarith.result = '@flatfieldedlist'
   saltred.saltarith.outpref=''
   saltred.saltarith.clobber = "yes"
   saltred.saltarith()

# Remove Cosmic Rays from Arcs:
saltred.saltcrclean(images='Reduced/*Arc*.fits',outpref='Reduced/cr.',crtype='fast',flux_ratio=0.25)

os.system('rm *.txt.*')
os.system('mv *.txt Reduced/')
os.system('chmod a+x calibrate.py')
os.system('mv calibrate.py Reduced/')
os.system('chmod a+x reduction.py')
os.system('mv reduction.py Reduced/')

#print >> reduction, "saltred.saltcrclean(images='tf.*.fits', outimages='', outpref='c', crtype='fast', thresh=5.0, mbox=5, bthresh=3.0, flux_ratio=0.1, bbox=11, gain=1.0, rdnoise=5.0, fthresh=5.0, bfactor=2, gbox=3, maxiter=5, multithread='no', clobber='no', logfile='salt.log', verbose='yes')"
#  Parameters suggested in: https://sciencewiki.salt.ac.za/index.php/SALT_Longslit_Reduction_Tutorial
print >> reduction, "saltred.saltcrclean(images='tf.*.fits', outpref='c', crtype='edge', bthresh=5.0, thresh=10.0, maxiter=5, mbox=25, flux_ratio=0.2, bbox=25, gain=1.0, rdnoise=5.0, fthresh=5.0, bfactor=2, gbox=3, multithread='no', clobber='yes', logfile='salt.log', verbose='yes')"

# Create residual images for inspection of Cosmic Rays removed
for j in range(len(targetfiles)):
   print >> reduction, "saltred.saltarith(operand1='tf.spectra"+str(j+1)+".fits', op='-', operand2='ctf.spectra"+str(j+1)+".fits', result='cr.residual"+str(j+1)+".fits', mode='h')"

# BACKGROUND subtraction
print >> reduction, "saltspec.specsky(images='ctf.*.fits', outimages='', outpref='b', method='normal', section='["+skybox+"]', clobber='yes', logfile='salt.log', verbose='yes')"

# Create 1D spectra in ASCII & FITS format from 2D spectra
for k in range(len(targetfiles)): 
   print >> reduction, "saltspec.specextract(images='bctf.spectra"+str(k+1)+".fits', outfile='bctf.spectra"+str(k+1)+".txt', method='normal', section='["+targetbox+"]', thresh=3.0, minsize=3.0, outformat='ascii', convert='yes', clobber='yes', logfile='salt.log', verbose='yes')"
   print >> reduction, "saltspec.specextract(images='bctf.spectra"+str(k+1)+".fits', outfile='bctf.spectra"+str(k+1)+".ms.fits', method='normal', section='["+targetbox+"]', thresh=3.0, minsize=3.0, outformat='FITS', convert='yes', clobber='yes', logfile='salt.log', verbose='yes')"
   # Use 2D spectra Without vacuum-correction or background-subtraction (to check the necessity of the air2vac correction using skylines)
   print >> reduction, "saltspec.specextract(images='ctf.spectra"+str(k+1)+".fits', outfile='ctf.spectra"+str(k+1)+".txt', method='normal', section='["+targetbox+"]', thresh=3.0, minsize=3.0, outformat='ascii', convert='yes', clobber='yes', logfile='salt.log', verbose='yes')"
   print >> reduction, "saltspec.specextract(images='ctf.spectra"+str(k+1)+".fits', outfile='ctf.spectra"+str(k+1)+".ms.fits', method='normal', section='["+targetbox+"]', thresh=3.0, minsize=3.0, outformat='FITS', convert='yes', clobber='yes', logfile='salt.log', verbose='yes')"
for kk in range(arccount):
   # Extract the 1D observed Arc
   print >> reduction, "saltspec.specextract(images='tcr.Arc"+str(kk+1)+".fits', outfile='tcr.Arc"+str(kk+1)+".txt', method='normal', section='["+targetbox+"]', thresh=3.0, minsize=3.0, outformat='ascii', convert='yes', clobber='yes', logfile='salt.log', verbose='yes')"
   print >> reduction, "saltspec.specextract(images='tcr.Arc"+str(kk+1)+".fits', outfile='tcr.Arc"+str(kk+1)+".ms.fits', method='normal', section='["+targetbox+"]', thresh=3.0, minsize=3.0, outformat='FITS', convert='yes', clobber='yes', logfile='salt.log', verbose='yes')"


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
print >> plots, "plot [] [0:] 'bctf.spectra"+str(1)+".txt' with lines title '"+timeobs+"    ("+targetfiles[0]+")'"
print >> plots, "unset label"
print >> plots, "unset title"
print >> plots, "set tmargin 0"
for l in range(1,len(targetfiles)-1):
   timeobs = str(imutil.hselect(images = targetfiles[l]+'[0]', fields = 'TIME-OBS', expr = 'yes',Stdout=1)).strip('[').strip(']').replace("'",'')
   print >> plots, "plot [] [0:] 'bctf.spectra"+str(l+1)+".txt' with lines title '"+timeobs+"    ("+targetfiles[l]+")'"
print >> plots, "set xlabel 'Wavelength [Angstrom]'"
print >> plots, "set xtics 0,200,10000"
if len(targetfiles) == 2:
   l = 0
if len(targetfiles) >= 2:
   timeobs = str(imutil.hselect(images = targetfiles[l+1]+'[0]', fields = 'TIME-OBS', expr = 'yes',Stdout=1)).strip('[').strip(']').replace("'",'')
   print >> plots, "plot [] [0:] 'bctf.spectra"+str(len(targetfiles))+".txt' with lines title '"+timeobs+"    ("+targetfiles[l+1]+")'"
print >> plots, "unset multiplot"
print >> plots, "exit"
plots.close()


######################
#  FIX -nan SPECTRA: #
######################

print >> reduction, 'import numpy'
print >> reduction, "os.system('ls bctf.spectra*txt > fixlist')"
print >> reduction, "filenames = numpy.loadtxt('fixlist',dtype='str')"

print >> reduction, "for i in range(len(filenames)):"
print >> reduction, "   fixedfile = open('fixedfile','w')"
print >> reduction, "   wavelength, flux, fluxerr = numpy.loadtxt(filenames[i],dtype='str', unpack='True')"
print >> reduction, "   for j in range(len(flux)):"
print >> reduction, "      if flux[j].count('-nan') == 0:"
print >> reduction, "         print >> fixedfile, wavelength[j], flux[j], fluxerr[j]"
print >> reduction, "      else:"
print >> reduction, "         print >> fixedfile, wavelength[j], 0, 0"
print >> reduction, "   fixedfile.close()"
print >> reduction, "   os.system('mv fixedfile '+filenames[i])"

print >> reduction, "os.system('gnuplot<plots')"
