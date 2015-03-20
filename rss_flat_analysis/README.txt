flatfld_rssimage.py lampflat skyflat fitslist (saveflat)

Flatfield correction for RSS, imaging mode

    lampflat: fits file name for daytime flat ("none" if not supplied)
    skyflat: fits file name for sky flat ("none" if not supplied)
    fitslist: unix commandline list of fits image files to be flattened.  A flattened 
        fits image file is written for each input file, with "f" prepended to the name. 
    option: if optional last argument is "saveflat", a fits file of the actual flat 
        field is saved for each image file.  "_flat" is appended to the file name.

You will need to edit the hardwired directory path "datadir" at line 16 of 
flatfld_rssimage.py. Then put imgdist.txt in datadir/spectrograph and 
track_1503_telflat.fits in datadir/telescope. 

There are three possible inputs to the flatfielding operation, the lamp and sky 
input flats, obtained for the observation, and the telescope track-flat model, 
supplied here as a fits hypercube which is grid of flats over track space. The 
lamp flat is used to correct high-frequency structure in the RSS detector, the 
sky flat is used to correct the RSS-dependent low-frequency structure (i.e., 
detector QE and gain), plus any time-dependent variations in the telescope 
vignetting not represented in the model), and the track grid models the 
low-frequency spatial dependence of the telescope vignetting, with its 
variation over track. The first two are assumed to rotate with the 
instrument, while the third does not. The flat process depends on what is 
supplied, lamp or sky, both, or neither.  In order of increasing preference,

1) Neither ("none none"). Only the telescope vignetting model is used, interpolated 
to the actual track position of the image. 
2) Lamp only. The telescope model in (1) is multiplied by the high frequency variations 
in the lamp flat. The lamp variations are evaluated separately for each CCD amplifier, 
so any difference in amplifier gains remains in the flat, and will be removed in the 
flattening.
3) Sky only. The sky image is processed to remove astronomical objects, and smoothed.
The smoothing is done by amplifier, so, again, the amplifier gains remain in the flat, 
and will be removed in flattening. This is then corrected by the ratio of telescope 
flat at the image track position to that at the track position of the sky flat.
4) Lamp and sky. The processing in (3) is done, and multiplied by the high frequency 
variations in the lamp flat. 

I have found that including a lamp flat is not very useful unless it is a combination
of multiple frames to increase the S/N.  The sky flat should not have any bright 
astronomical objects (especially galaxies), and is best taken with the clear filter
to maximize sky S/N and minimize filter ghosts (a 20 second exposure is plenty).

This is a "beta" version of this software. In particular, the track model (version
"1503") only goes out to 1.6 m from the center of track (not into the "corners"), 
and is assumed to be circularly symmetric in track space, which is not quite true. A 
better model, requiring more sky data, is planned. 

Below is sample output from the routine, with no lamp flat, a sky flat from the 
beginning of a track, applied to 7 images in a track, with the optional save of 
the flats themselves. This takes about 30 seconds on my machine, with most of 
the time spent at the beginning setting up the track model. 

-khn-oderbolz-55-> flatfld_rssimage.py none cmbxgpP201304280014.fits c*1[5-9].fits c*2[0-1].fits saveflat
RSS Imaging Flatfield using telescope model  1503
saved : fcmbxgpP201304280015.fits cmbxgpP201304280015_flat.fits
saved : fcmbxgpP201304280016.fits cmbxgpP201304280016_flat.fits
saved : fcmbxgpP201304280017.fits cmbxgpP201304280017_flat.fits
saved : fcmbxgpP201304280018.fits cmbxgpP201304280018_flat.fits
saved : fcmbxgpP201304280019.fits cmbxgpP201304280019_flat.fits
saved : fcmbxgpP201304280020.fits cmbxgpP201304280020_flat.fits
saved : fcmbxgpP201304280021.fits cmbxgpP201304280021_flat.fits

