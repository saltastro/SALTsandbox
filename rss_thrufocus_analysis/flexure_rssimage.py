#! /usr/bin/env python

# Compute RSS flexure vs rho, any imaging mask

import os, sys, time
import numpy as np
from datetime import datetime
from astropy.io import fits as pyfits

import warnings
import matplotlib
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
matplotlib.use('PDF')
plt.ioff()

from pyraf import iraf
from iraf import pysalt
from saltobslog import obslog
from thrufoc_rssimage import sextract

np.set_printoptions(threshold=np.nan)

datadir = "/d/carol/Synched/software/SALT/pipeline/"    

# ---------------------------------------------------------------------------------
def flexure_rssimage(fitslist,option=""):

    print str(datetime.now()), "\n"
    if option == "filesave": prefix = raw_input("file prefix: ")

    pixel = 15.                                                     # pixel size in microns
    pix_scale=0.125 
    sexparams = ["X_IMAGE","Y_IMAGE","FLUX_ISO","FLUX_MAX","FLAGS","CLASS_STAR",    \
                "X2WIN_IMAGE","Y2WIN_IMAGE","XYWIN_IMAGE","ERRX2WIN_IMAGE"]
    np.savetxt("qred_thrufoc.param",sexparams,fmt="%s")
    fmaxcol,flagcol,xvarcol,yvarcol,xerrcol = (3,4,6,7,9)           # column nos (from 0) of data in sextractor
    imagestooclosefactor = 3.0                                      # too close if factor*sep < sqrt(var)
    gaptooclose = 1.25                                              # arcsec
    edgetooclose = 1.25                                             # arcsec
    rattolerance = 0.25    
    toofaint = 250.                                                 # FMAX counts
    galaxydelta = 0.4                                               # arcsec
    MOSimagelimit = 1.                                              # arcsec
    deblend = .005                                                  # default

    flexposns = len(fitslist)
    obsdict=obslog(fitslist)

    image_f = [fitslist[fpos].split(".")[0][-12:] for fpos in range(flexposns)]
    dateobs = int(image_f[0][:8])
    if dateobs > 20110928:  rho_f = np.array(obsdict["TRKRHO"]).astype(float)
    else:                   rho_f = np.array(obsdict["TELRHO"]).astype(float)
    catpos = np.argmin(np.abs(rho_f))

    cbin,rbin = np.array(obsdict["CCDSUM"][catpos].split(" ")).astype(int)
    maskid =  obsdict["MASKID"][catpos].strip() 
    if maskid == 'P000000N99': colspots = 21
    else: colspots = 1
    filter =  obsdict["FILTER"][catpos].strip()
    grating =  obsdict["GRATING"][catpos].strip()
    rows,cols = pyfits.getdata(fitslist[catpos]).shape
    isspec = (obsdict["GR-STATE"][catpos][1] =="4")

    print "\nMask:      ", maskid
    print "Filter:    ", filter
    print "Grating:   ", grating

    if isspec:
        print "Use flexure_rssspec for spectral flexure"
        exit()

#   make catalog of good spots using image closest to rho=0 (capital S after cull)

    sex_js = sextract(fitslist[catpos],deblend=deblend)
             
    fluxmedian = np.median(np.sort(sex_js[2])[-10:])    # median of 10 brightest                               
    ok_s = (sex_js[2] > fluxmedian/100.)                # cull bogus spots 
    sexcols = sex_js.shape[0]
    Spots = ok_s.sum()
    sexdata_jfS = np.zeros((sexcols,flexposns,Spots))
    sexdata_jfS[:,catpos] = sex_js[:,ok_s]
    xcenter = 0.5*(sexdata_jfS[0,catpos].min() +  sexdata_jfS[0,catpos].max())
    ycenter = 0.5*(sexdata_jfS[1,catpos].min() +  sexdata_jfS[1,catpos].max())  

    print "\n     fits      rho  spots  rshift  cshift  rslope  cslope  rmserr "
    print   "               deg         arcsec  arcsec  arcmin  arcmin   bins"
    print ("%12s %5.1f %5i "+5*"%7.2f ") % \
        (image_f[catpos], rho_f[catpos], Spots, 0., 0., 0., 0., 0.)

    if option == "filesave":
        np.savetxt(prefix+"Spots.txt",sexdata_jfS[:,catpos].T,   \
            fmt=2*"%9.2f "+"%9.0f "+"%9.1f "+"%4i "+"%6.2f "+3*"%7.2f "+"%11.3e")

#   find spots in flexure series, in order of increasing abs(rho), and store sextractor output

    row_fd = np.zeros((flexposns,2))
    col_fd = np.zeros((flexposns,2))

    for dirn in (1,-1):
        refpos = catpos
        posdirlist = np.argsort(dirn*rho_f)
        poslist = posdirlist[dirn*rho_f[posdirlist] > rho_f[refpos]]

        for fpos in poslist:
            col_S,row_S = sexdata_jfS[0:2,refpos,:]   
            sex_js = sextract(fitslist[fpos],"sexwt.fits",deblend=deblend)     
            bintol = 16/cbin                                  # 2 arcsec tolerance for finding star
            binsqerr_sS = (sex_js[1,:,None] - row_S[None,:])**2 + (sex_js[0,:,None] - col_S[None,:])**2
            S_s = np.argmin(binsqerr_sS,axis=1)
        # First compute image shift by averaging small errors
            rowerr_s = sex_js[1] - row_S[S_s]
            colerr_s = sex_js[0] - col_S[S_s]
            hist_r,bin_r = np.histogram(rowerr_s,bins=32,range=(-2*bintol,2*bintol))    
            drow = rowerr_s[(rowerr_s > bin_r[np.argmax(hist_r)]-bintol) & \
                (rowerr_s < bin_r[np.argmax(hist_r)]+bintol)].mean()
            hist_c,bin_c = np.histogram(colerr_s,bins=32,range=(-2*bintol,2*bintol))    
            dcol = colerr_s[(colerr_s > bin_c[np.argmax(hist_c)]-bintol) & \
                (colerr_s < bin_c[np.argmax(hist_c)]+bintol)].mean()
        # Now refind the closest ID
            binsqerr_sS = (sex_js[1,:,None] - row_S[None,:] -drow)**2 + \
                (sex_js[0,:,None] - col_S[None,:] -dcol)**2
            binsqerr_s = binsqerr_sS.min(axis=1)
            isfound_s = binsqerr_s < bintol**2
            S_s = np.argmin(binsqerr_sS,axis=1)
            isfound_s &= (binsqerr_s == binsqerr_sS[:,S_s].min(axis=0))
            isfound_S = np.array([S in S_s[isfound_s] for S in range(Spots)])
            sexdata_jfS[:,fpos,S_s[isfound_s]] = sex_js[:,isfound_s]
            drow_S = sexdata_jfS[1,fpos]-sexdata_jfS[1,catpos]
            dcol_S = sexdata_jfS[0,fpos]-sexdata_jfS[0,catpos]
            row_fd[fpos],rowchi,d,d,d = np.polyfit(sexdata_jfS[0,catpos,isfound_S]-xcenter,   \
                 drow_S[isfound_S],deg=1,full=True)
            col_fd[fpos],colchi,d,d,d = np.polyfit(sexdata_jfS[1,catpos,isfound_S]-ycenter,   \
                 dcol_S[isfound_S],deg=1,full=True)
            rms = np.sqrt((rowchi+colchi)/(2*isfound_S.sum()))

            print ("%12s %5.0f %5i "+5*"%7.2f ") % (image_f[fpos], rho_f[fpos], isfound_S.sum(), \
                row_fd[fpos,1]*rbin*pix_scale, col_fd[fpos,1]*cbin*pix_scale, \
                60.*np.degrees(row_fd[fpos,0]),-60.*np.degrees(col_fd[fpos,0]), rms)
            if option == "filesave":
                np.savetxt(prefix+"flex_"+str(fpos)+".txt",np.vstack((isfound_S,drow_S,dcol_S)).T,  \
                    fmt = "%2i %8.3f %8.3f")
                np.savetxt(prefix+"sextr_"+str(fpos)+".txt",sexdata_jfS[:,fpos].T)

#   make plots

    fig,plot_s = plt.subplots(2,1,sharex=True)

    plt.xlabel('Rho (deg)')
    plt.xlim(-120,120)
    plt.xticks(range(-120,120,30)) 
    fig.set_size_inches((8.5,11))
    fig.subplots_adjust(left=0.175)

    plot_s[0].set_title(str(dateobs)+" Imaging Flexure") 
    plot_s[0].set_ylabel('Mean Position (arcsec)')
    plot_s[0].set_ylim(-0.5,4.)
    plot_s[1].set_ylabel('Rotation (arcmin ccw)')      
    plot_s[1].set_ylim(-10.,6.)

    plot_s[0].plot(rho_f,row_fd[:,1]*rbin*pix_scale,marker='D',label='row')
    plot_s[0].plot(rho_f,col_fd[:,1]*rbin*pix_scale,marker='D',label='col')
    plot_s[1].plot(rho_f,60.*np.degrees(row_fd[:,0]),marker='D',label='row')
    plot_s[1].plot(rho_f,-60.*np.degrees(col_fd[:,0]),marker='D',label='col')

    plot_s[0].legend(fontsize='medium',loc='upper center')                     
    plotfile = str(dateobs)+'_imflex.pdf'
    plt.savefig(plotfile,orientation='portrait')

    if os.name=='posix':
        if os.popen('ps -C evince -f').read().count(plotfile)==0: os.system('evince '+plotfile+' &')

    os.remove("out.txt")
    os.remove("qred_thrufoc.param") 
    os.remove("sexwt.fits")                           
    return

# ---------------------------------------------------------------------------------

if __name__=='__main__':
    fitslist=sys.argv[1:]
    option = ""
    if (fitslist[-1] == "filesave"):
        option = fitslist.pop()
    flexure_rssimage(fitslist,option)

# debug:
# cd /d/pfis/khn/20110804/sci
# flexure_rssimage.py m*E3??.fits
# cd /d/pfis/khn/20160725/sci
# flexure_rssimage.py m*0059.fits m*006[0-7].fits filesave
