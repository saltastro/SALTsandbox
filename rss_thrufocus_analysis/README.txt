thrufoc_rsscartesian.py fitslist (fwhmfile= debug=)

Thru-focus analysis for a cartesian mask (eg N99 or N02).

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
thrufoc_rssspec.py. Then put qred.sex there, create a subdirectory "spectrograph" 
there, and put spec.txt, gratings.txt, and specfocus.txt in that subdirectory. 

The samples 20141122_1800Ne_thrufoc.txt and 20141122_2300Ne_thrufoc.txt show you what I 
get when I run it on longslit thru-focus runs taken on 20141122 with the Ne lamp 
through the 1800 and 2300 l/mm gratings. On my machine, this takes less than one 
minute. The program looks at the center image in the focus run and finds the brightest 
lines in each of the six ccd amplifiers, derives a wavelength for each using the 
spectrograph model (good enough to be within a few Angstroms), and predicts the best 
focus for that line using the grating, wavelength, and temperature, and the same focus 
model used on PCON. Then it looks at each focus image and derives the line width as a 
function of column, printing out the mean width for each line, so you can see the rough 
best-focus. Each line now has a focus curve as a function of column, which is 
block-averaged and spline-interpolated into a smooth curve. The best focus, tip and tip 
curvature of the focal curve of each line is computed from this, and, comparing those, 
the tilt and tilt curvature. The tilt is relative to the *expected* tilt derived from 
the spectrograph model, so all grating configurations should nominally give the same 
tilt. The exception to this is the 2300 l/mm grating, which has appreciable VPH-induced 
focal power, casing its spectral tilt to deviate from the other gratings by about 4 
arcmin. 

At the very bottom, advice is given on the screw adjustments for the detector which 
would zero out the tilt and tip errors.
