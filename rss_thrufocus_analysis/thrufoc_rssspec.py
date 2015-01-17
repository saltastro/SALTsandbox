#! /usr/bin/env python

# Compute RSS focal surface for six lines across spectrum in a spectral thrufocus

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

datadir = "/d/carol/Synched/software/SALT/pipeline/"    

# ---------------------------------------------------------------------------------
def blksmooth1d(ar,blk,mask=[]):
# blkaverage (ignoring zeros), then spline interpolate result

    blks = ar.shape[0]/blk
    offset = (ar.shape[0] - blks*blk)/2
    arr = np.copy(ar); arr[mask] = 0.
    arrr = arr[offset:offset+blks*blk]                  # cut to integral blocks
    sum = arrr.reshape(blks, blk).sum(axis=1)

    count = (arrr.reshape(blks, blk) > 0).sum(axis=1)
    arblk = np.zeros_like(sum)
    okcounts = (count > blk/2)
    arblk[okcounts] = sum[okcounts]/count[okcounts]
    gridblk = np.arange(offset+blk/2,offset+blks*blk+blk/2,blk)
    grid = np.arange(ar.shape[0])
    arsmooth = ip.griddata(gridblk[okcounts], arblk[okcounts], grid, method="cubic", fill_value=0.)

    return arsmooth

# ---------------------------------------------------------------------------------
def rsslam(grating, grang, artic, dpix_c):
#   compute on-axis wavelengths for array of columns, using full spectrograph model
#   grating: grating name
#   grang: grating angle (deg)
#   artic: articulation angle (deg)
#   dpix_c: array of column coordinates, in unbinned pixels relative to center

    spec=np.loadtxt(datadir+"/spectrograph/spec.txt",usecols=(1,))
    Grat0,Home0,ArtErr,T2Con,T3Con=spec[0:5]
    FCampoly=spec[5:11]
    grname=np.loadtxt(datadir+"/spectrograph/gratings.txt",dtype=str,usecols=(0,))
    grlmm,grgam0=np.loadtxt(datadir+"/spectrograph/gratings.txt",usecols=(1,2),unpack=True)
    grnum = np.where(grname==grating)[0][0]
    lmm = grlmm[grnum]
    alpha_r = np.radians(grang+Grat0)
    beta0_r = np.radians(artic*(1+ArtErr)+Home0)-alpha_r
    gam0_r = np.radians(grgam0[grnum])
    lam0 = 1e7*np.cos(gam0_r)*(np.sin(alpha_r) + np.sin(beta0_r))/lmm
    ww = lam0/1000. - 4.
    fcam = np.polyval(FCampoly,ww)
    disp = (1e7*np.cos(gam0_r)*np.cos(beta0_r)/lmm)/(fcam/.015)
    dfcam = 3.162*disp*np.polyval([FCampoly[x]*(5-x) for x in range(5)],ww)
    T2 = -0.25*(1e7*np.cos(gam0_r)*np.sin(beta0_r)/lmm)/(fcam/47.43)**2 + T2Con*disp*dfcam
    T3 = (-1./24.)*3162.*disp/(fcam/47.43)**2 + T3Con*disp
    T0 = lam0 + T2 
    T1 = 3162.*disp + 3*T3
    X = dpix_c/3162.
    lam_c = T0+T1*X+T2*(2*X**2-1)+T3*(4*X**3-3*X)

#   print artic, grang, grnum, lmm, lam0, ww, fcam, disp, dfcam
#   print T0, T1, T2, T3
    return lam_c

# ---------------------------------------------------------------------------------
def thrufoc_rssspec(fitslist,option=""):

    if option == "filesave": prefix = raw_input("file prefix: ")

    imgdist=np.loadtxt(datadir+"/spectrograph/imgdist.txt",usecols=(1,2))
    Fsclpoly=imgdist[4: ,0]
    wv_w, focoffs_w = np.loadtxt(datadir+"/spectrograph/specfocus.txt",unpack=True) # focus vs wavelength
    Fm = 736.-330.
    dFta, dFto = 3.5, -34                                                       # Focus vs temperature coefficients
    grating_g = np.array(["PG0300","PG0900","PG1300","PG1800","PG2300","PG3000"])
    dFg_g = np.array([0,0,0,0,135.,43.])                                        # grating fopcus offsets
    pixel = 15.                                                                 # pixel size in microns
    pixarcsec = pixel/Fsclpoly[5]
    turnfac_d = [8.0,6.7]                                                       # arcmin/turn for tip/tilt adustment
    turnmicr = 200.                                                             # microns/turn
    focposns = len(fitslist)

#   find 6 useful lines from center focus (left,right of each ccd)
    lines = 6
    col_l = np.zeros(lines+1)
    pix1_l = np.array([100.,900.,2200.,3000.,4300.,5100.])                      # unbinned pixels
    pix2_l = pix1_l + 800.
    catposn = focposns/2
    hdulist_image = pyfits.open(fitslist[catposn])
    rows,cols = hdulist_image["SCI"].data.shape
    row_r = np.arange(rows)
    hdr = hdulist_image[0].header
    rbin,cbin = map(int,(hdr["CCDSUM"].split(" ")))
    grating =  (hdr["GRATING"]).strip()
    grang = float(hdr["GRTILT"])
    artic =  float(hdr["CAMANG"])
    mask = (hdr["MASKID"]).strip()
    lamp = (hdr["LAMPID"]).strip()
    if "COLTEM" in hdr:                                                         # allow for annoying SALT fits version change
        coltem = float(hdr["COLTEM"]) 
        camtem = float(hdr["CAMTEM"])
    else:
        coltem = float(hdr["COLTEMP"]) 
        camtem = float(hdr["CAMTEMP"])     

    print "Grating:        ", grating
    print "Artic (deg):    ", artic
    print "Gr Angle (deg): ", grang
    print "Mask:           ", mask
    print "Lamp:           ", lamp
    print "Col Temp:       %7.1f" % (coltem)
    print "Cam Temp:       %7.1f" % (camtem)
    onaxis_c = np.copy(hdulist_image["SCI"].data)[rows/2-32/cbin:rows/2+32/cbin,:].mean(axis=0)
    for line in range(lines): 
        col_l[line] = np.argmax(onaxis_c[pix1_l[line]/cbin : pix2_l[line]/cbin]) + pix1_l[line]/cbin
    col_l[lines] = cols/2
    wav_l = rsslam(grating, grang, artic, (col_l - cols/2)*cbin)
    focoff_l = ip.interp1d(wv_w, focoffs_w,kind ='cubic')(wav_l)
    gratno = np.where(grating_g==grating)[0]
    predfoc_l = Fm + focoff_l + dFto*(coltem-20.) + dFta*(camtem-20.) + dFg_g[gratno]
    print "\n line         0      1      2      3      4      5      center"
    print " column  ", (lines*"%6.0f "+"%9.0f") % tuple(col_l)
    print " wavelen ", (lines*"%6.0f "+"%9.0f") % tuple(wav_l)
    print " predfoc ", (lines*"%6.0f "+"%9.0f") % tuple(predfoc_l), "\n"

#   calculate width(row) for each chosen line in each focus image
    print "\n image   focus  rows     ------------line width(arcsec)------------"
    print "        (micr)           0       1       2       3       4       5 " 
    image_f = np.zeros((focposns))
    focus_f = np.zeros((focposns))
    allgoodrow_r = np.ones((rows)).astype(np.bool_)
    width_flr = np.zeros((focposns,lines,rows))

    for fpos in range(focposns):   
        image_f[fpos] = int(fitslist[fpos].split(".")[0][-4:])
        hdulist_image = pyfits.open(fitslist[fpos])
        focus_f[fpos] = float(hdulist_image[0].header["FO-POS-V"])
        for line in range(lines):
            slice_rc = np.copy(hdulist_image["SCI"].data)[:,col_l[line]-100/cbin:col_l[line]+300/cbin]
            argmax_r = np.maximum(20/cbin,np.minimum(380/cbin,np.argmax(slice_rc,axis=1)))          # column position of line
            line_rc = ((np.arange(400/cbin)[None,:]>argmax_r[:,None]-20/cbin) \
                 & (np.arange(400/cbin)[None,:]<argmax_r[:,None]+20/cbin))
            slice_rc = slice_rc[line_rc].reshape(rows,(40/cbin-1))                  # line is straightened out
            slice_rc -= slice_rc.min(axis=1).mean()                                 # remove background
            flux_r = np.sum(slice_rc,axis=1)
            fmax_r = slice_rc.max(axis=1)
            fmax = np.sort(fmax_r)[-64]                                             # typical line strength 
            goodrow_r = (fmax_r > fmax/10.) & (fmax_r < 1.5*fmax) & (flux_r > 0.)   # get rid of spikes and ends of line                                                     
            allgoodrow_r &= goodrow_r
            width_flr[fpos,line,goodrow_r] = flux_r[goodrow_r]/fmax_r[goodrow_r]

        print "%5i %7.0f %5i" % (image_f[fpos], focus_f[fpos], allgoodrow_r.sum()),
        width_lr = width_flr[fpos]
        width_l =  width_lr[:,allgoodrow_r].mean(axis=1)
        print lines*"%7.2f " % tuple(width_l*cbin*pixarcsec)

#   compute best focus(row) for each line

    goodrow_lr = ((width_flr > 0.).sum(axis=0) == focposns)
    focusgrid_F = np.arange(focus_f[0],focus_f[-1],1.)
    bestfocussmooth_lr = np.zeros((lines,rows))
    tip_l = np.zeros(lines)
    tipcurv_l = np.zeros(lines)

    print "\n line    focus     tip   tip curvature\n         (micr)  (arcmin)   (micr)"

    for line in range(lines):
        bestfocus_r = np.zeros(rows)
        focuscurve_r = np.zeros(rows)
        goodrow_r = ((width_flr[:,line,:] > 0.).sum(axis=0) == focposns)
        for r in range(rows):
            if goodrow_r[r]: 
                width_f = width_flr[:,line,r]; 
                focuscurve_F = ip.interp1d(focus_f,width_flr[:,line,r],kind ='cubic')(focusgrid_F) 
                bestfocus_r[r] = focusgrid_F[np.argmin(focuscurve_F)]
        meanfoc = bestfocus_r[goodrow_r].mean()
        goodrow_r &= (np.absolute(bestfocus_r - meanfoc) < 150.)
        bestfocus_r[np.logical_not(goodrow_r)] = 0.
        tip_l[line],meanfoc= np.polyfit(row_r[goodrow_r],bestfocus_r[goodrow_r],1)
        bestfocussmooth_lr[line,:] = blksmooth1d(bestfocus_r,64)
        goodrow_r = (bestfocussmooth_lr[line,:] > 0.)
        focuscurve_r[goodrow_r] = bestfocussmooth_lr[line,goodrow_r] - np.polyval((tip_l[line],meanfoc),row_r[goodrow_r])
        center_r = goodrow_r & (row_r> 0.4*rows) & (row_r< 0.6*rows)
        edge_r = goodrow_r & ((row_r< 0.15*rows) | (row_r> 0.85*rows))
        tipcurv_l[line] = focuscurve_r[center_r].mean() - focuscurve_r[edge_r].mean()

        print "%6.0f %7.0f  %7.2f  %7.0f" % (wav_l[line],meanfoc,60*np.degrees(tip_l[line]/(pixel*rbin)),tipcurv_l[line])

#   compute focal plane tilt relative to predicted
    tilt,dfoc0= np.polyfit(col_l[0:lines],(bestfocussmooth_lr[:,rows/2]-predfoc_l[0:lines]),1)
    dfoc0 += tilt*cols/2
    tiltcurv = bestfocussmooth_lr[2:4,rows/2].mean() - bestfocussmooth_lr[0:6:5,rows/2].mean()

    print "\n centerfoc  slope (arcmin)  curvature (micr)"
    print "   (micr)    tip     tilt     tip   tilt "
    print "%8.1f %8.2f %8.2f %6.0f %6.0f" % (predfoc_l[lines]+dfoc0,60*np.degrees(tip_l.mean()/(pixel*rbin)),  \
                                        60*np.degrees(tilt/(pixel*rbin)),tipcurv_l[2:4].mean(),tiltcurv)
    print "\n Suggested adjustment turns"
    print " tip:   top :", "%5.1f" % (60*np.degrees(tip_l.mean()/(pixel*rbin))/turnfac_d[0])
    print "     bottom :", "%5.1f" % (-60*np.degrees(tip_l.mean()/(pixel*rbin))/turnfac_d[0])
    print " tilt: left :", "%5.1f" % (-60*np.degrees(tilt/(pixel*rbin))/turnfac_d[1])
    print " mean:  all :", "%5.1f" % (dfoc0/turnmicr)

    if option=="filesave":
        focplanefile = open(prefix+"_focplane.txt","w")
        print >> focplanefile, " row   focshift------------->"
        for r in range(rows): print >> focplanefile, ("%5i " % (r)), (lines*"%8.2f " % tuple(bestfocussmooth_lr[:,r]))

    return

# ---------------------------------------------------------------------------------
if __name__=='__main__':
    fitslist=sys.argv[1:]
    option = ""
    if (fitslist[-1] == "filesave"):
        option = fitslist.pop()

    thrufoc_rssspec(fitslist,option)

