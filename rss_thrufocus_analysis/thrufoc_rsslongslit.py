#! /usr/bin/env python

# Compute RSS focal tip and rotation for longslit imaging thrufocus

import os, sys, time
import numpy as np
from scipy import interpolate as ip
import pyfits
from datetime import datetime
from thrufoc_rssspec import blksmooth1d

np.set_printoptions(threshold=np.nan)

datadir = "/d/carol/Synched/software/SALT/pipeline/" 

# ---------------------------------------------------------------------------------
def thrufoc_rsslongslit(fitslist,option=""):
    """Analyse imaging through-focus data of longslit 

    Parameters
    ----------
    fitslist: list
        List of filenames 

    option: str
        filesave:  save text file of rotation and focus tilt

    """
    """
    _f: index of focus run
    _r: image row (bins)
    _c: image column (bins)
    _F: index of culled focus run

    """

    print str(datetime.now()),'\n'
    if option == "filesave": prefix = raw_input("file prefix: ")

    pixel = 15.                                             # pixel size in microns
    pixarcsec = 0.125                                       # arcsec/pix
    turnfac_d = [8.0,6.7]                                   # arcmin/turn for tip/tilt adustment
    turnmicr = 200.                                         # microns/turn
    focposns = len(fitslist)

#   find slit at middle of focus run
    catposn = focposns/2
    hdulist_image = pyfits.open(fitslist[catposn])
    hdr = hdulist_image[0].header
    cbin,rbin = map(int,(hdr["CCDSUM"].split(" ")))
    rows,cols = hdulist_image["SCI"].data.shape
    mask = (hdr["MASKID"]).strip()
    lamp = (hdr["LAMPID"]).strip()
    filter =  (hdr["FILTER"]).strip()
    if "COLTEM" in hdr:                            # allow for annoying SALT fits version change
        coltem = float(hdr["COLTEM"]) 
        camtem = float(hdr["CAMTEM"])
    else:
        coltem = float(hdr["COLTEMP"]) 
        camtem = float(hdr["CAMTEMP"])     
    slitarcsecs = float(mask[2:5])/10.
    slitbins = slitarcsecs*8./cbin
    onaxis_c = hdulist_image["SCI"].data[rows/2-32/rbin:rows/2+32/rbin,:].mean(axis=0)
    col = np.median(onaxis_c.argsort()[-(2+slitbins):])
    row_r = np.arange(rows)

#   predict best focus
    wav_w, focoffs_w = np.loadtxt(datadir+"/spectrograph/specfocus.txt",unpack=True) 
    filter_f = np.loadtxt(datadir+"/spectrograph/filters.txt",dtype='str',usecols=(0,))
    wav_f,dFf_f = np.loadtxt(datadir+"/spectrograph/filters.txt",usecols=(1,2),unpack=True)
    fno = np.argwhere(filter==filter_f)[0]
    focoff = ip.interp1d(wav_w, focoffs_w,kind ='cubic')(wav_f[fno]) + dFf_f[fno]
    Fm = 736.-330.
    dFta, dFto = 3.5, -34  
    predfoc = Fm + focoff + dFto*(coltem-20.) + dFta*(camtem-20.)

    print "Mask:   ", mask, "  slitwidth: ", slitarcsecs, slitbins, " arcsec, bins"
    print "Filter:    ", filter
    print "Lamp:   ", lamp
    print "Col Temp:       %7.1f" % (coltem)
    print "Cam Temp:       %7.1f" % (camtem)
    print "Predicted Foc:  %7.1f" % (predfoc)

#   calculate slit width vs row in each focus image
    print "\n image  focus  rows line width(arcsec)"
    image_f = np.zeros((focposns),dtype=int)
    focus_f = np.zeros((focposns))
    allgoodrow_r = np.ones((rows)).astype(np.bool_)
    width_fr = np.zeros((focposns,rows))
    col_fr = np.zeros((focposns,rows))

    for fpos in range(focposns):   
        image_f[fpos] = int(fitslist[fpos].split(".")[0][-4:])
        hdulist_image = pyfits.open(fitslist[fpos])
        focus_f[fpos] = float(hdulist_image[0].header["FO-POS-V"])
        slice_rc = np.copy(hdulist_image["SCI"].data)[:,col-100/cbin:col+100/cbin]
        slice_rc -= slice_rc.min(axis=1).mean()                                # remove background
        col_fr[fpos] = np.median(slice_rc.argsort(axis=1)[:,-(2+slitbins):],axis=1)   # column position vs row
        flux_r = np.sum(slice_rc,axis=1)
        fmax_r = slice_rc.max(axis=1)
        fmax = np.sort(fmax_r)[-64]                                             # typical line strength 
        goodrow_r = (fmax_r > fmax/10.) & (fmax_r < 1.5*fmax) & (flux_r > 0.)   # get rid of spikes and ends of line                                                     
        allgoodrow_r &= goodrow_r
        width_fr[fpos,goodrow_r] = flux_r[goodrow_r]/fmax_r[goodrow_r]

        print "%5i %7.0f %5i" % (image_f[fpos], focus_f[fpos], allgoodrow_r.sum()),
        width_r = width_fr[fpos]
        width =  width_r[allgoodrow_r].mean()
        print "%7.2f " % (width*cbin*pixarcsec)

#   compute best focus as a function of row
    goodrow_r = ((width_fr > 0.).sum(axis=0) == focposns)
    focusgrid_F = np.arange(focus_f[0],focus_f[-1],1.)
    bestfocussmooth_r = np.zeros(rows)

    bestfocus_r = np.zeros(rows)
    focuscurve_r = np.zeros(rows)
    goodrow_r = ((width_fr[:,:] > 0.).sum(axis=0) == focposns)
    for r in range(rows):
        if goodrow_r[r]: 
            width_f = width_fr[:,r]
            focuscurve_F = ip.interp1d(focus_f,width_fr[:,r],kind ='cubic')(focusgrid_F) 
            bestfocus_r[r] = focusgrid_F[np.argmin(focuscurve_F)]
    meanfoc = bestfocus_r[goodrow_r].mean()
    goodrow_r &= (np.absolute(bestfocus_r - meanfoc) < 150.)
    bestfocus_r[np.logical_not(goodrow_r)] = 0.
    bestfpos = np.argmin(np.abs(meanfoc-focus_f))

#   fit slit position vs row for best image, best focus vs row to fin slope and errors
    (rot,d),((vrot,d),(d,d)) = np.polyfit(row_r[goodrow_r],col_fr[bestfpos,goodrow_r],1,cov=True)
    (tip,meanfoc),((vtip,d),(d,d))= np.polyfit(row_r[goodrow_r],bestfocus_r[goodrow_r],1,cov=True)
    if option == "filesave":
        np.savetxt(prefix+"_colfocus.txt",np.vstack((goodrow_r,col_fr[bestfpos],bestfocus_r)).T, \
            fmt="%2i %4i %6.0f",header="goodrow col bestfocus") 
    bestfocussmooth_r = blksmooth1d(bestfocus_r,64)
    goodrow_r = (bestfocussmooth_r > 0.)
    focuscurve_r[goodrow_r] = bestfocussmooth_r[goodrow_r] - np.polyval((tip,meanfoc),row_r[goodrow_r])
    center_r = goodrow_r & (row_r> 0.4*rows) & (row_r< 0.6*rows)
    edge_r = goodrow_r & ((row_r< 0.15*rows) | (row_r> 0.85*rows))
    tipcurv = focuscurve_r[center_r].mean() - focuscurve_r[edge_r].mean()

    print "\n best image: ",image_f[bestfpos]
    print "\n focus    rot     tip  tip curvature\n (micr) (arcmin) (arcmin) (micr)"
    print "%6.0f  %6.2f  %6.2f  %6.0f" \
        % (meanfoc,60*np.degrees(rot*(cbin/rbin)),60*np.degrees(tip/(pixel*rbin)),tipcurv)
    print " +/-    %6.2f  %6.2f" \
        % (60*np.degrees(np.sqrt(vrot)*(cbin/rbin)),60*np.degrees(np.sqrt(vtip)/(pixel*rbin)))

    print "\n Suggested adjustment turns"
    print " tip:   top :", "%5.1f" % (60*np.degrees(tip/(pixel*rbin))/turnfac_d[0])
    print "     bottom :", "%5.1f" % (-60*np.degrees(tip/(pixel*rbin))/turnfac_d[0])

    return

# ---------------------------------------------------------------------------------
if __name__=='__main__':
    fitslist=sys.argv[1:]
    option = ""
    if (fitslist[-1] == "filesave"):
        option = fitslist.pop()
    thrufoc_rsslongslit(fitslist,option)
