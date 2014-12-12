#! /usr/bin/env python

# Compute RSS imaging flatfield for a specific image

import os, sys, time
import numpy as np
from scipy.optimize import fmin
from scipy.ndimage import interpolation as nd
from scipy import interpolate as ip
import pyfits
from math import pi

np.set_printoptions(threshold=np.nan)

class Blackhole:
    """Utility class used to suppress some output
    """
    def write(self, *msg):
        pass
STDOUT = sys.stdout

# ---------------------------------------------------------------------------------
def sextract(fits,option=""):
# run sextractor to find objects
    global hdulist_image, shape_d
    hdulist_image = pyfits.open(fits)
    shape_d = np.asarray(hdulist_image["SCI"].data.shape)
    rcbin_d = np.array(hdulist_image[0].header["CCDSUM"].split(" ")).astype(int)
    image_rc = np.copy(hdulist_image["SCI"].data)
    sigma = 5.                    
    pix_scale=0.125
    r_ap=2.0/(pix_scale*rcbin_d.min())         
    sat=0.99*image_rc.max()
    cmd='sex %s -c $PYTHONPATH/qred.sex -CATALOG_NAME %s -DETECT_THRESH %f -PHOT_APERTURES %f -SATUR_LEVEL %f' \
            % (fits, "out.txt", sigma, r_ap, sat)
    os.system(cmd+" &> /dev/null")
    sex_js = np.loadtxt("out.txt").transpose()
    os.remove("out.txt")
    return sex_js

# ---------------------------------------------------------------------------------
def thrufoc_rssimage(name,fitslist,option=""):
    global sex_js, rd, B, pixarcsec, gridsize, niter
    global focus_f, rc0_dgg, usespots_gg, fwhminterp_g                                            # for cube analysis

    focposns = len(fitslist)

#   make catalog of grid spots using middle image
    pixel = 15.                                                                  # pixel size in microns
    gridsize = 21
    turnfac_d = [8.0,6.7]                                                        # arcmin/turn for tip/tilt adustment
    rcoff_df = np.zeros((2,focposns))    
    catposn = focposns/2
    rcbin_d = np.array(pyfits.getheader(fitslist[catposn],0)["CCDSUM"].split(" ")).astype(int)
    maskid =  (pyfits.getheader(fitslist[catposn],0)["MASKID"]).strip() 
    focuscat = pyfits.getheader(fitslist[catposn],0)["FO-POS-V"] 
    sex_js = sextract(fitslist[catposn],option)
    sexcols = sex_js.shape[0]

#   get first guess positions from distortion model
    datadir = "/d/carol/Synched/software/SALT/pipeline/"    
    gridsep = 22.                                                               # arcsec spacing of mask grid
    fittol = 3e-3
    imgdist=np.loadtxt(datadir+"/spectrograph/imgdist.txt",usecols=(1,2))

    rcd_d = imgdist[1,::-1]                                                     # note backwards from file
    rd, A,B = imgdist[2,0],imgdist[3,0],imgdist[3,1]
    Fsclpoly=imgdist[4: ,0]
    pixarcsec = pixel/Fsclpoly[5]

    if      maskid == "P000000N99":  rc0_d = (-56.,+12.)
    elif    maskid == "P000000N01":  rc0_d = (-22.,+18.)
    else:
        print maskid, " is not a thrufocus mask"
        exit()

#   optimize model
    x0 = np.zeros((7))    
    x0[0:2] = rc0_d 
    x0[2:4] = rcd_d 
    x0[4:7] = A,0.,gridsep

    print "Distortion Fit for focus ", focuscat
    print "    pix grid0x   grid0y   dist0x   dist0y      A   rot(deg) gridsep(arcsec) " 
    print (" start"+4*"%8.2f "+2*"%8.4f "+"%8.2f") % tuple(x0) 

    niter = 0
    rms = poserr(x0)
    rcoff_d = ((sex_js[(1,0),s_gg[gridsize/2,gridsize/2]] - rc_dgg[:,gridsize/2,gridsize/2]))     # center up the center spot first
    x0[0:2] += rcoff_d                  
    sys.stdout = Blackhole()                                    # suppress output of the fmin function
    xopt = fmin(poserr, x0, xtol=fittol)
    sys.stdout = STDOUT                                         # Restore normal output
    rc0_dgg = np.copy(rc_dgg)                                   # save predicted grid positions for fitted model

    print (" fit  "+4*"%8.2f "+2*"%8.4f "+"%8.2f") % tuple(xopt)  

    image_f = np.zeros((focposns))                                                                   
    focus_f = np.zeros((focposns)) 
    sexdata_jfgg = np.zeros((sexcols,focposns,gridsize,gridsize))
    foundcount_gg = np.zeros((gridsize,gridsize))

#   using fitted distortion model, find grid spots in focus series, and store sextractor output
    print "\n image focus  spots"
    for fpos in range(focposns):   
        sex_js = sextract(fitslist[fpos])
        image_f[fpos] = int(fitslist[fpos].split(".")[0][-4:])
        rms = poserr(xopt)
        rcoff_d = ((sex_js[(1,0),s_gg[gridsize/2,gridsize/2]] - rc_dgg[:,gridsize/2,gridsize/2]))    # center up the center spot first
        x = xopt; x[0:2] = xopt[0:2] + rcoff_d
        rms = poserr(x)
        focus_f[fpos] = pyfits.getheader(fitslist[fpos],0)["FO-POS-V"]
        sexdata_jfgg[:,fpos,found_dgg[0],found_dgg[1]] = sex_js[:,s_gg[found_dgg]]
        print "%5i %6.0f %5i" % (image_f[fpos], focus_f[fpos], found_dgg[0].shape[0])

    foundcount_gg = (sexdata_jfgg[0] > 0).sum(axis=0)
    sexdata_jfgg *= (foundcount_gg > (focposns - 3))[None,None,:,:]                # throw out cr's, etc
    foundcount_gg = (sexdata_jfgg[0] > 0).sum(axis=0)

#    for g in range(gridsize): print "%3i " % (g), gridsize*"%3i " % tuple(foundcount_gg[g,:])

#   fit a plane to best-fwhm to get best focus and focal plane tilt 
#   all fully-populated grid points for imaging, central cross of grid for longslit
#   start with focus for lowest fwhm of central grid point 
    x0 = np.array([focus_f[np.argmin(sexdata_jfgg[10,:,gridsize/2,gridsize/2])],0.,0.])
    args_gg = np.indices((gridsize,gridsize))

    print "\nBest focal plane          arcmin        turns in adj "
    print   "         focus(micr)    tip    tilt      top   right  fwhm(arcsec) spots"

    niter = 0
    usespots_gg = (foundcount_gg == focposns) & ((args_gg[0] == gridsize/2) | (args_gg[1] == gridsize/2))
    fwhmdata_fg = sexdata_jfgg[10][:,usespots_gg]
    fwhminterp_g = np.empty(usespots_gg.sum(),dtype='object')
    for g in range(usespots_gg.sum()): fwhminterp_g[g] = ip.interp1d(focus_f,fwhmdata_fg[:,g],kind ='cubic')
    sys.stdout = Blackhole()                                                       # suppress output of the fmin function
    xls = fmin(fwhmnet, x0, xtol=fittol)
    rms = pixarcsec*fwhmnet(xls)
    arcmin_d = 60.*np.degrees(xls[1:]/pixel)
    sys.stdout = STDOUT  
    print " longslit %8.1f " % xls[0], 2*"%7.2f " % tuple(arcmin_d), \
            2*"%7.2f " % tuple(arcmin_d/turnfac_d), "%8.2f %8i" % (rms, usespots_gg.sum())

    niter = 0
    usespots_gg = (foundcount_gg == focposns)
    fwhmdata_fg = sexdata_jfgg[10][:,usespots_gg]
    fwhminterp_g = np.empty(usespots_gg.sum(),dtype='object')
    for g in range(usespots_gg.sum()): fwhminterp_g[g] = ip.interp1d(focus_f,fwhmdata_fg[:,g],kind ='cubic')
    sys.stdout = Blackhole()                                                       
    ximg = fmin(fwhmnet, x0, xtol=fittol)
    rms = pixarcsec*fwhmnet(ximg)
    arcmin_d = 60.*np.degrees(ximg[1:]/pixel)
    sys.stdout = STDOUT                                                            # Restore normal output     
    print " imaging  %8.1f " % ximg[0], 2*"%7.2f " % tuple(arcmin_d), \
            2*"%7.2f " % tuple(arcmin_d/turnfac_d), "%8.2f %8i" % (rms, usespots_gg.sum())

#   write out fits cubes for selected data
    colout_i = (8,9,10)
    outfile_i = ("theta","ellip","fwhm")

    print "\n Saving ",outfile_i, " data cubes"

    for i in range(len(colout_i)):
        outfile = name+"_"+outfile_i[i]+".fits"
        hdulist_image["SCI"].data = sexdata_jfgg[colout_i[i]]
        hdulist_image.writeto(outfile, clobber = "True")                        
    return

# ---------------------------------------------------------------------------------
def poserr(x):
    global rms, rc_dgg, s_gg, found_dgg, rcoff_d, niter

    rc0_d = x[0:2]
    rcd_d = x[2:4]
    A,rot,gridsep = x[4:7]

#   predict grid positions based on current model

    pixtol = 0.5*gridsep/pixarcsec                                               # pixel tolerance for finding spot

    ea_dgg = gridsep*(np.indices((gridsize,gridsize)) - gridsize/2)/pixarcsec    # undistorted grid in pixels
    c = np.cos(np.radians(rot))
    s = np.sin(np.radians(rot))
    rotate = np.transpose([[c, s],[-1.*s, c]])

    rc_dgg = np.transpose(np.dot(np.transpose(ea_dgg,(1,2,0)),rotate),(2,0,1)) + rc0_d[:,None,None]
    r_gg = np.sqrt((rc_dgg[0]-rcd_d[0])**2 + (rc_dgg[1]-rcd_d[1])**2)
    rc_dgg = rcd_d[:,None,None] + (rc_dgg - rcd_d[:,None,None])*(1. + A*(r_gg[None,:,:]/rd)**2 + B*(r_gg[None,:,:]/rd)**4) \
        + shape_d[:,None,None]/2.

# identify sextractor spots and compute position err    
    s_gg = np.argmin((sex_js[1,:,None,None] - rc_dgg[0,None,:,:])**2 + (sex_js[0,:,None,None] - rc_dgg[1,None,:,:])**2,axis=0)
    found_dgg = np.where((sex_js[1,s_gg] - rc_dgg[0,:,:])**2 + (sex_js[0,s_gg] - rc_dgg[1,:,:])**2 < pixtol**2)
    sumsq = ((sex_js[1,s_gg[found_dgg]] - rc_dgg[0][found_dgg])**2 + (sex_js[0,s_gg[found_dgg]] - rc_dgg[1][found_dgg])**2).sum()
    rms = np.sqrt(sumsq)/found_dgg[0].shape[0]

#    print ("%4i %8.3f %4i " % (niter, rms, found_dgg[0].shape[0])), (x.shape[0]*"%11.4f " % tuple(x))
 
    niter += 1
    return rms
# ---------------------------------------------------------------------------------
def fwhmnet(x):
    global niter, fplane_g, fwhm_g

    fplane_g = (x[0] + x[1]*(rc0_dgg[0]- rc0_dgg[0,gridsize/2,gridsize/2]) \
                     + x[2]*(rc0_dgg[1]- rc0_dgg[1,gridsize/2,gridsize/2]))[usespots_gg]                            # location of focal plane
    fplane_g = np.maximum(fplane_g,focus_f[0])
    fplane_g = np.minimum(fplane_g,focus_f[-1])
    fwhm_g = np.zeros_like(fplane_g)
    for g in range(usespots_gg.sum()): fwhm_g[g] = fwhminterp_g[g](fplane_g[g])
    rms = np.sqrt((fwhm_g**2).sum()/usespots_gg.sum())                                                              # rms net fwhm
#    print niter, "%8.3f " % (rms), ("%8.2f "+2*"%8.6f ") % tuple(x)
    niter += 1
    return rms

# ---------------------------------------------------------------------------------
if __name__=='__main__':

    name=sys.argv[1]
    fitslist=sys.argv[2:]
    option = ""
    if (fitslist[-1] == "debug"):
        option = fitslist.pop()

    thrufoc_rssimage(name,fitslist,option)
