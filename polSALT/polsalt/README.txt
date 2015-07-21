reducepoldata.py obsdate
========================

Reduce RSS spectropolarimetric data through image, wavecal, extract, and stokes 
    evaluation steps.

Copy the script 'reducepoldata.py' to the directory that has the data that you want to 
reduce. Edit this file to point to the polsalt directory, and the proposal code. 
Comment out steps you don't need (eg, have already performed).

The directory where you will run the data should have the format such that in the top 
level directory there should be the obsdate (e.g. 20061030) and then the raw directory 
inside that. The raw directory should only include the files that you want reduced so 
any extra files which are sometimes included in your directory should be removed. 

Calibration data for this command are stored in the polsalt/data directory. The command 
does the following steps:

imred.  Basic image reduction.
------------------------------ 
Logged in 'im'+obsdate+'.log'. Modeled after imred in zsalt (with variance and bad 
pixel layers, required for polarimetry), except as noted: 

prepare. (same as zsalt) 
bias. (same as zsalt)
variance. Uses master unbinned bad pixel mask stored in /data , and bins it according 
   to fits header. 
gain. (same as zsalt) 
crosstalk. (same as zsalt) 
flatfield. Not performed. 
cosmic ray removal. (same as zsalt, skipped for lamp data) 
mosaic. Fixed errors in mosaicing bpm and variance. 
cleanup. (same as zsalt, leaves m*.fits file) 

Basic polarimetry steps below are logged in 'specpol'+obsdate+'.log'. Non-polarimetric 
images (based on beamsplitter out) are skipped.  For now, the assumption is made that all
images are of the same object/configuration (where "configuration" is defined by 
grating, grating angle, articulation, and waveplate pattern). This will be fixed in the 
future. 
 
specpolmap. Split beams, wavelength calibrate.
----------------------------------------------
The arc image is split into O and E arcs, and the beamsplitter chromaticism is 
removed to straighten the dispersion temporarily. Each is fed to specidentify in 
interactive mode. The result is used to make an image which is a map of wavelength for 
the original (unstraightened) O and E beams. Each image in the input list (with its 
VAR and BPM extensions) is now split into O and E images, which are stacked in the 
third fits dimension. The wavelength map is added to each as a "WAV" extension. A prefix
"w" is put on the image name. There are two stokes layers, O and E, each with SCI, VAR, 
BPM, and WAV extensions. 

specpolextract. Define target extraction aperture, and extract spectra.
-----------------------------------------------------------------------
All non-arc images are combined (O and E separately).  An extraction aperture is 
determined by finding the cross-dispersion maximum and going out to where the signal has 
dropped to 2% of maximum. A background aperture then runs from this point to the end of 
the slit, or to where the signal rises back to 2% of max (this assumes we are looking 
at the central object in a long slit). The O and E spectra are then extracted (and VAR 
and BPM updated) using these apertures for each image. The same extraction aperture is 
used for all images in order to avoid systematic noise between the different 
observations in the stokes pattern. The resulting spectra are stored with a "e" prefix. 
There are two stokes layers, O and E, each with SCI, VAR, BPM, and WAV extensions. The y 
axis in the fits files is virtual (length 1). 

specpolrawstokes. Combine O,E spectra in filter pairs to produce "raw stokes" spectra.
--------------------------------------------------------------------------------------
Makes a table locating each "filter pair" in a polarization pattern (eg, contiguous 
exposures with halfwave positions differing by 45 degrees, which flip E and O). 
Use WAV map to resample E and O onto a uniform wavelength grid, so they may be combined 
to make an intensity and a raw stokes polarization spectrum. Save the results as files 
with name "object_config_hnm_r.fits", where n and m are the halfwave stations and r is 
the repeat number, like QTH1-QTH2_c0_h04_01.fits. For now, only linear polarization is 
supported. When circular and all-stokes is added, these will append a "qnm" to the "hnm" 
to defined the filterpair. There are two stokes layers, I and S(unnormalized), each 
with SCI, VAR,and BPM extensions. 

specpolfinalstokes.  Combine raw stokes spectra and apply calibration for final result.
---------------------------------------------------------------------------------------
(For future: put in track-dependent corrections to raw stokes spectra here). All raw 
stokes spectra for a particular object/config are combined, first combining repeats of a 
particular filter pair, then, depending on the waveplate pattern, combining different 
pairs. For Linear-HI, this allows calculation of a "chisq" assessment of any errors 
beyond readout and photon noise.  They should be near 1.  If very different from that, 
some data culling may be required. Now apply polarimetric efficiency and instrumental 
position angle calibrations, and, for sky observations, add in the rotator and telescope 
azimuth PA correction. Save as "object_config_stokes.fits" file. For linear patterns, 
there are three stokes layers, I, Q,and U (unnormalized), each with SCI, VAR,and BPM 
extensions. There is a fourth VAR layer for the QU covariance, which is needed to compute 
P and PA errors for rotated linear polarization. For Circular patterns there will be 2 
layers, I and V, and for All-Stokes, 4 layers, I, Q, U, V, with a fifth VAR layer for 
QU covariance. The calibration file is identified in a "POLCAL" FITS keyword. It is in 
"polcal.txt" in the /data directory, and is identified in its second line. 

______________________________________________________________________________________

specpolview fitslist binoption saveoption.  Plot and text analysis of stokes files.
===================================================================================
    fitslist: one or more raw or final stokes files (but not both).
    binoption: either not given (unbinned), nnb (use bins of width nn bins), or nn% (bin 
        to nn %S or %P error, for raw or linear polarization spectra)
    saveoption: either "plot" (save plot as pdf), "text" (save text table of what is in 
        plot), or both (eg "plottext")

Plot intensity and stokes spectra vs wavelength in stacked plots.  Polarization is in % 
(normalized), and PA is in degrees. If there is more than one fits, they are overplotted, 
with a legend identifying which is which. If there is binning, the bins are plotted as 
error symbols, x being the bin width and y being the stokes error. 


