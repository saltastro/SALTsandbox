Descriptions for:
thrufoc_rsslongslit.py: Imaging thrufocus of longslit
thrufoc_rsscartesian.py: Imaging thrufocus of cartesian mask (eg N99 or N01)
thrufoc_rssimage.py: Imaging thrufocus of any imaging mask (or no mask)
thrufoc_rssspec.py: Longslit spectral thrufocus
flexure_rss.py: Imaging flexure run
--------------------------------------------------------------------------------------

thrufoc_rsslongslit.py fitslist (filesave)

Imaging thru-focus analysis of longslit 

    fitslist: unix commandline list of fits files
    option : if optional last argument = "filesave", a text file "(prefix)_colfocus.txt"
      is saved containing a "good row" flag, the column position, and the best focus of 
      the slit as a function of row.  The file name prefix is requested. 

You will need to edit hardwired directory path "datadir" at line 14 of 
thrufoc_rsslongslit.py. Then put qred.sex there.

The sample 20111113_longslit_thrufoc.txt shows you what I get when I run it on the 0.6 
arcsec longslit thru-focus data taken on 20111113 with the QTH lamp. On my machine, this 
takes 5-10 secs. The configuration and the predicted best focus ("Autofocus" on PCON) are
printed out.  The program looks at the center image in the focus run to find the slit. 
Then it looks at each focus image and derives the slit width as a function of row, 
printing out the mean width for each line, so you can see the rough best-focus. The slit 
now has a focus curve as a function of row, which is block-averaged and spline-
interpolated into a smooth curve. The best focus, and the tip and tip curvature of the 
focal plane is computed from this. For the image closest to best focus, the column 
position of the slit as a function of row is fit with a line to obtain the rotation of the 
slit relative to the vertical.  The best focus, rotation, tip and tip curvature of the 
slit (with errors) is then printed out. 

At the very bottom, advice is given on the screw adjustments for the detector which 
would zero out the tip error.
--------------------------------------------------------------------------------------

thrufoc_rsscartesian.py fitslist (fwhmfile= debug=)

Thru-focus analysis for a cartesian mask (eg N99 or N01).

    fitslist: unix commandline list of fits files 
    optional: fwhmfile (string). Saves fwhmfile_fwhm.fits cube of fwhm vs row, col,focus
    optional: debug=True.  Saves sextractor file for center of focus run.    

You will need to edit hardwired directory path "datadir" at line 23 of 
thrufoc_rsscartesian.py. Then put qred.sex there.

The sample 20141212_6645_thrufoc_rsscartesian.txt shows you what I 
get when I run it on the 20141212 6645 filter thrufocus data taken with the "N99" 
cartesian mask. On my machine, this takes about 30 seconds. 

The slitmask is identified, and the known tip and tilt of that mask is given, which is 
removed from the detector plane tip, tilt to be calculated. The filter and temperature 
data are for information only, and are not used by the program. The program first looks at 
the middle image in the run (presumably near focus) using SExtractor, identifies the 
spots with a cartesian grid using a first-guess distortion model ("start"), and computes a 
best-fit distortion model ("fit"). It then runs SExtractor on all the images to tabulate 
image properties of the spots as a function of position and focus. "spots" is the total 
number of spots found in each image. A focus curve is computed for each spot.  A best 
focus position for each focus curve is computed by minimizing the rms difference to the 
median focus curve as a function of focal shift. Finally, it fits the best focus as 
a function of position to a polynomial surface to obtain the mean central focus, tip, and 
tilt, and the mean fwhm for three cases:  First, for longslit applications it fits a plane 
to spots along the x and y axes.  Second, for imaging it fits a plane to all spots. Third 
("maskrim"), it fits all spots to a second-order polynomial surface to obtain a mask 
curvature and the tilt and tip of the focal plane at the edge of the field of view (ie, 
the orientation of the mask frame).

At the very bottom, advice is given on the screw adjustments for the detector which 
would zero out the tilt and tip errors.
 
--------------------------------------------------------------------------------------
thrufoc_rssimage.py fitslist (filesave)

Thru-focus analysis for any imaging mask (or no mask)

    fitslist: unix commandline list of fits files 
    option : if optional last argument = "filesave", text files are saved: 
      "(prefix)_cull.txt" contains the cull status of each star in each image (see source 
      for key), and "(prefix)_focplane.txt" contains the focus shift for each good focus 
      curve. The file name prefix is requested. 

You will need to edit hardwired directory path "datadir" at line 14 of 
thrufoc_rssimage.py. Then put qred.sex there. 

The samples 20141124_6645_N99_thrufoc.txt and 20111105_open_thrufoc.txt show you what I 
get when I run it on the 20141124 6645 filter thrufocus data taken with the "N99" 
cartesian mask and on 20111105 thrufocus on sky (no mask) (Note that the 1105 run has 
an incorrect indication of mask of PL0060N001 - the fits header is incorrect and should 
read P00000N00. The mask, filter, and temperature data at the top are for information 
only, and are not used by the program). On my machine, this takes about 4 minutes. 
Basically, the program looks at all the images in the run using SExtractor, then 
identifies those spots with the spots found in the middle image to make a map of focus 
curves. "spots" is the total number of spots found in each image. "rshift" and "cshift" 
is the row and column shift in arcsec of each image from the central image (note the 
guiding drift in the open mask run). "xfwhm", "yfwhm", and "fwhm" is the median fwhm in 
arcsec for each image, so you can see the rough position of best focus. "rmserr" is the 
rms error in bins in the spot positions, relative to the central image. The spots are 
then put though a culling process to remove close doubles, objects too faint, 
odd-shaped objects, ones too near the gaps and the edge of field, and finally, objects 
which have fewer than half the focus points in a curve or which are abnormally big 
(i.e., galaxies in the open mask thru-focus). Then it gets a focus shift for each focus 
curve by minimizing the rms difference to the median focus curve as a function of focal 
shift. Finally, it fits the focus shift as a function of position with a plane to 
obtain the mean central focus, tip, and tilt (with errors), then with a second-order 
polynomial to obtain a curvature and the tilt and tip of the focal plane at the edge of 
the field of view (ie, the orientation of the mask frame). 

No advice on detector alignment is now given, since I find that the focal plane 
orientation for imaging with masks is dominated by the mask mounting error.

--------------------------------------------------------------------------------------
thrufoc_rssspec.py fitslist (filesave)

Thru-focus analysis for longslit spectroscopy

    fitslist: unix commandline list of fits files
    option : if optional last argument = "filesave", a text file "(prefix)_focplane.txt"
      is saved containing the smoothed focus shift for each line.  The file name prefix 
      is requested. 

    You will need to edit hardwired directory path "datadir" at line 22 of 
thrufoc_rssspec.py to be a directory of your choice. Then put qred.sex there, create a 
subdirectory "spectrograph" there, and put the spec_yyymmdd.txt, gratings.txt, 
imgdist.txt, and specfocus.txt in that subdirectory. 

The samples 20141122_1800Ne_thrufoc.txt and 20141122_2300Ne_thrufoc.txt show you what I 
get when I run it on longslit thru-focus runs taken on 20141122 with the Ne lamp 
through the 1800 and 2300 l/mm gratings. On my machine, this takes less than one 
minute. The program looks at the center image in the focus run and finds the brightest 
lines in each of the six ccd amplifiers, derives a wavelength for each using the 
spectrograph model (good enough to be within a few Angstroms), and predicts the best 
focus for that line using the grating, wavelength, and temperature, and the same focus 
model used on PCON. Then it looks at each focus image and derives the line width as a 
function of row, printing out the mean width for each line, so you can see the rough 
best-focus. Each line now has a focus curve as a function of row, which is block-averaged 
and spline-interpolated into a smooth curve. The best focus, tip and tip curvature of the 
focal curve of each line is computed from this, and, comparing those, the tilt and tilt 
curvature. The tilt is relative to the *expected* tilt derived from the spectrograph 
model, so all grating configurations should nominally give the same tilt. The exception to 
this is the 2300 l/mm grating, which has appreciable VPH-induced focal power, casing its 
spectral tilt to deviate from the other gratings by about 4 arcmin. 

At the very bottom, advice is given on the screw adjustments for the detector which 
would zero out the tilt and tip errors.

--------------------------------------------------------------------------------------
flexure_rssimage.py fitslist (filesave)

Imaging flexure analysis 

    fitslist: unix commandline list of fits files
    option : if optional last argument = "filesave", text files "(prefix)_flex_(fpos).txt"
      and "(prefix)_sextr_(fpos).txt" are saved for each flexure image (0,...,fpos).  The 
      "flex" files containing a "good row" flag, and the row and column of each object 
      relative to the rho=0 image (in pixels).  The sextr files contain the detailed 
      SeXtract output for each image. The file name prefix is requested. 

    You will need to edit hardwired directory path "datadir" at line 25 of 
flexure_rssimage.py to be a directory of your choice. Then put qred.sex there. 

The sample 20160725_imflex.txt shows you what I get when I run it on the cartesian mask 
(N99) imaging flexure data taken on 20160725 with the QTH lamp. On my machine, this takes 
just over a minute.  The program looks at the image with rho closest to 0 and locates the 
mask spots (attempting to cull out flaws and cosmic rays), and uses this as a baseline.  
Then it proceeds out from rho = 0, finding where the spots go, while limiting motion to 2 
arcsec to avoid spoilers. For each image it determines the mean shift in row and column 
of the whole mask image, and the rotation of the mask image, using the row and column 
motion separately to estimate changing distortion.  Using these 4 parameters, it 
determines the remaining residual rms in bins, which assesses the remaining distortion 
change.  Finally, it generates a pdf graphic of the four parameters, saving them in the 
.pdf files. 

--------------------------------------------------------------------------------------
flexure_rssspec.py imagefits specfitslist (filesave)

Spectral flexure analysis 

    imagefits: fits image of the mask (rho=0 image from imaging flexure run is good)
    specfitslist: unix commandline list of spectral flexure fits files
    option : if optional last argument = "filesave", text files "(prefix)_flex_(fpos).txt"
      and "(prefix)_sextr_(fpos).txt" are saved for each flexure image (0,...,fpos).  The 
      "flex" files containing a "good row" flag, and the row and column of each object 
      relative to the rho=0 image (in pixels).  The sextr files contain the detailed 
      SeXtract output for each image. The file name prefix is requested.

    You will need to edit hardwired directory path "datadir" at line 29 of 
flexure_rssspec.py to be a directory of your choice. Then put qred.sex there, create a 
subdirectory "spectrograph" there, and put the spec_yyymmdd.txt files, and gratings.txt
in that subdirectory. 

The sample 20160725_grflex.txt shows you what I get when I run it on the cartesian mask 
(N99) spectral flexure data taken on 20160725 with the QTH lamp. On my machine, this takes 
just over a minute.  The program first maps the mask, reporting the total number of holes,
the number of mask rows and columns, and the pixel position of the center.  Then it looks
at the spectral flexure image closest to rho = 0 to identify the spectral lines, each of 
which generates an image which looks like the mask.  This is culled down to the brightest 
lines that form relatively complete mask images, listed in a table.  The overall shift 
caused by the grating at rho=0 is shown, together with the line identification rms, in 
pixels.  The identified spots are then identified in each image in the flexure series, 
working out from rho = 0 as in the imaging flexure analysis above. The four flexure 
parameters are then evaluated for each line, the mean row and column shift and the 
apparent rotation in row and column about the center of the mask for that line, and these 
are printed out, with the rms of the fit to this 4-parameter model.  Finally, it generates 
a pdf graphic of the four parameters as a function of rho and line, saving them in a .pdf 
file.  Column and row are shown together, as diamond and square markers.
