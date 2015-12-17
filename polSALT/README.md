polSALT
=======

Spectropolarimetric reductions for SALT longslit observations.  


Instructions
============

1. Clone the git repository onto your device

   git clone https://github.com/saltastro/SALTsandbox/polSALT.git

2. Copy the script 'reducepoldata.py' to the directory that has the
   data that you want to reduce.   Edit this file to point to the 
   polsalt directory.

3. The directory where you will run the data should have the format
   such that in the top level directory there should be the
   obsdate (e.g. 20061030) and then the raw directory inside that.
   The raw directory should only include the files that you want
   reduced so any extra files which are sometimes included in your
   directory should be removed.

4. Run reducepoldata.py with

   python reducepoldata.py 20061030

   Replace the obsdate with the appropriate observation date

5. Reducepoldata should step through all the tasks needed to produce
   wavelength calibrated stokes spectra.   There will be two interactive
   steps for each arc (O and E beams) to identify arc lines.

6. Once complete, inspect the final product in the sci directory
   that was created.  There should be "raw stokes" files 
   "name_config_hnm.fits" for each halfwaveplate pair nm in the pattern,
   plus a "name_config_stokes.fits" for the calibrated stokes spectra.
   (see README.txt in polsalt for details). The stokes spectra are
   in the third dimension of the fits file, so DS9 will bring up
   a cube slice widget to select which is displayed.  For "raw stokes"
   there are two slices, I and "S" (the stokes polarization corresponding 
   to the filter pair), and for linear polarization calibrated stokes (the
   only pattern currently supported), three slices, I, Q, and U.
   These are unnormalized stokes parameters.

7. To examine stokes spectra, use specpolview.py. You will either need the 
   full path to specpolview.py or copy it into the directory. Examples:  

   To overplot the intensity and normalized raw stokes I, S/I for 
   objectname, configuration c0:
   python specpolview.py (objectname)_c0_h*.fits

   To plot I, P/I, PA (deg) for this object/configuration:
   python specpolview.py (objectname)_c0_stokes.fits

   Polarization spectra may be binned by wavelength and to constant error, 
   and the plot and/or an ascii table of what is in the plot may be saved
   using options described in polsalt/README.txt.



