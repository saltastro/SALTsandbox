#! /usr/bin/env python

# Compute RSS spectral flexure vs rho

import os, sys, time
import numpy as np
from datetime import datetime
from scipy.interpolate import interp1d
from scipy import interpolate as ip
from scipy import linalg as la
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
from thrufoc_rssspec import rsslam, datedfile
np.set_printoptions(threshold=np.nan)
np.seterr(invalid='raise')

datadir = "/d/carol/Synched/software/SALT/pipeline/"    

# ---------------------------------------------------------------------------------
def flexure_rssspec(imagefits,fitslist,option=""):

    print str(datetime.now())

    if option == "filesave": prefix = raw_input("\nFile prefix: ")
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

    imagehdr = pyfits.getheader(imagefits)
    if imagehdr["GR-STATE"][1] == "4":
        print "First fits file "+imagefits+" is not image of mask"
        exit()

    flexposns = len(fitslist)
    obsdict=obslog(fitslist)

    image_f = [fitslist[fpos].split(".")[0][-12:] for fpos in range(flexposns)]
    dateobs = obsdict["DATE-OBS"][0].replace("-","")
    if int(dateobs) > 20110928: rho_f = np.array(obsdict["TRKRHO"]).astype(float)
    else:                       rho_f = np.array(obsdict["TELRHO"]).astype(float)
    
    catpos = np.argmin(np.abs(rho_f))
    cbin,rbin = np.array(obsdict["CCDSUM"][catpos].split(" ")).astype(int)
    maskid =  obsdict["MASKID"][catpos].strip() 
    filter =  obsdict["FILTER"][catpos].strip()
    grating =  obsdict["GRATING"][catpos].strip()
    rows,cols = pyfits.getdata(fitslist[catpos]).shape
    isspec = (obsdict["GR-STATE"][catpos][1] =="4")
    if not isspec:
        print "Use flexure_rssimage for image flexure analysis"
        exit()
    grang = float(obsdict["GRTILT"][catpos])
    artic = float(obsdict["CAMANG"][catpos])
    lamp = obsdict["LAMPID"][catpos].strip()

    print "\nMask:           ", maskid
    print "Filter:         ", filter
    print "Grating:        ", grating
    print "Artic (deg):    ", artic
    print "Gr Angle (deg): ", grang
    print "Lamp:           ", lamp

#   map the mask spots _m using the imaging fits file
    sex_js = sextract(imagefits,deblend=deblend)
    flux_s = sex_js[2]
    fluxmedian = np.median(np.sort(flux_s)[-10:])
    okm_s = (flux_s > fluxmedian/10)               # cull bogus spots
    maskholes = okm_s.sum()
    r_m = sex_js[1,okm_s]
    c_m = sex_js[0,okm_s]

#   find mask rows _R, tabulate 
    histr_b, binr_b = np.histogram(r_m,bins=rows/10,range=(0,rows))
    bin0_R = np.where((histr_b[1:]>0) & (histr_b[:-1]==0))[0]
    bin1_R = np.where((histr_b[1:]==0) & (histr_b[:-1]>0))[0]
    maskRows = bin0_R.shape[0]
    bin_m = np.digitize(r_m,binr_b) - 1
    R_m = np.array([np.where((bin_m[m] >= bin0_R) & (bin_m[m] <= bin1_R))[0][0] \
                    for m in range(maskholes)])

#   find mask cols _C, tabulate 
    histc_b, binc_b = np.histogram(c_m,bins=cols/10,range=(0,cols))
    bin0_C = np.where((histc_b[1:]>0) & (histc_b[:-1]==0))[0]
    bin1_C = np.where((histc_b[1:]==0) & (histc_b[:-1]>0))[0]
    maskCols = bin0_C.shape[0]
    bin_m = np.digitize(c_m,binc_b) - 1
    C_m = np.array([np.where((bin_m[m] >= bin0_C) & (bin_m[m] <= bin1_C))[0][0] \
                    for m in range(maskholes)])

#   identify mask center = optical axis
    if maskid == 'P000000N99':              # symmetric mask
        Raxis = maskRows/2
        Caxis = maskCols/2
    elif maskid == 'P000000N03':            # mask with centered cross 
        Raxis = np.where((np.argmax(histr_b) >= bin0_R) & (np.argmax(histr_b) <= bin1_R))[0][0]  
        Caxis = np.where((np.argmax(histc_b) >= bin0_C) & (np.argmax(histc_b) <= bin1_C))[0][0]
    else:
        print "Not a valid flexure mask"
        exit()
    maxis = np.where((R_m==Raxis)&(C_m==Caxis))[0][0]
    raxis = r_m[maxis]
    caxis = c_m[maxis]

    print "\nMask_Holes Rows Cols  r axis   c axis \n                      pixels   pixels"
    print "  %5i %5i %5i %8.1f %8.1f" % (maskholes,maskRows,maskCols,raxis*rbin,caxis*cbin)

#    np.savetxt(dateobs+'_'+"mask.txt",np.vstack((r_m,c_m,sex_js[2,okm_s],R_m)).T,fmt="%10.2f")

#   get linelist, predict spots in spectral image
    wavcent = rsslam(grating, grang, artic, 0.,dateobs)
    specfile = datedfile(datadir+"spectrograph/spec_yyyymmdd.txt",dateobs)
    FCampoly=np.loadtxt(specfile,usecols=(1,))[5:11]
    fcam = np.polyval(FCampoly,(wavcent/1000. - 4.))
    lampfile=iraf.osfn("pysalt$data/linelists/"+lamp+".salt")
    wav_l,int_l = np.loadtxt(lampfile,unpack=True)

    maxdalpha = -np.degrees((cols/2)*cbin*pixel/(1000.*fcam))
    maxgamma = np.degrees((rows/2)*rbin*pixel/(1000.*fcam))
    maxwav = rsslam(grating,grang,artic, cols*cbin/2,dateobs,-maxdalpha,0)
    minwav = rsslam(grating,grang,artic,-cols*cbin/2,dateobs, maxdalpha,maxgamma)
    ok_l = (wav_l >= minwav) & (wav_l <= maxwav)
    wav_l = wav_l[ok_l]
    int_l = int_l[ok_l]
    lines = wav_l.shape[0]
    col_ml = np.zeros((maskholes,lines))
    dcol_c = np.arange(-(cols*cbin/2),(cols*cbin/2))
    for m in range(maskholes):
        dalpha = -np.degrees((c_m[m]-caxis)*cbin*pixel/(1000.*fcam))
        gamma = np.degrees((r_m[m]-raxis)*rbin*pixel/(1000.*fcam))
        wav0,wav1 = rsslam(grating,grang,artic,dcol_c[[0,-1]],dateobs,dalpha,gamma=gamma)
        ok_l = ((wav_l > wav0) & (wav_l < wav1))
        colwav = interp1d(rsslam(grating,grang,artic,dcol_c,    \
                dateobs,dalpha=dalpha,gamma=gamma), dcol_c)
        col_ml[m,ok_l] = colwav(wav_l[ok_l]) + caxis*cbin 

#    np.savetxt(dateobs+"_col_ml.txt",np.vstack((R_m,C_m,col_ml.T)),fmt="%8.1f")     

#   identify mask hole and wavelength for spots in spec image closest to rho=0
    os.remove("sexwt.fits")   
    sex_js = sextract(fitslist[catpos],"",deblend=deblend)
    r_s = sex_js[1]
    c_s = sex_js[0]
    flux_s = sex_js[2]
    spots = r_s.shape[0]
    fluxmedian = np.median(np.sort(sex_js[2])[-10:])
    ok_s = (flux_s > fluxmedian/30)                    # cull bogus spots
             
#   find spectral bin rows RR in candidates R0, cull non-spectra
    histr_b, binr_b = np.histogram(r_s[ok_s],bins=rows/10,range=(0,rows))
    histr_b[[0,-1]] = 0
    bin0_R0 = np.where((histr_b[1:]>0) & (histr_b[:-1]==0))[0] + 1
    bin1_R0 = np.where((histr_b[1:]==0) & (histr_b[:-1]>0))[0]
    bin_s = np.digitize(r_s,binr_b) - 1

    maxcount_R0 = np.array([(histr_b[bin0_R0[R0]:bin1_R0[R0]+1]).max() \
                        for R0 in range(bin0_R0.shape[0])])
    ok_R0 = (maxcount_R0 > 3)                              
    specrows = ok_R0.sum()                              # cull down to spectra RR
    bin0_RR = bin0_R0[ok_R0]
    bin1_RR = bin1_R0[ok_R0]
    ok_s &= ((bin_s >= bin0_RR[:,None]) & (bin_s <= bin1_RR[:,None])).any(axis=0)
    RR_s = -np.ones(spots)
    r_RR = np.zeros(specrows)
    for RR in range(specrows):
        isRR_s = ok_s & np.in1d(bin_s,np.arange(bin0_RR[RR],bin1_RR[RR]+1))
        RR_s[isRR_s] = RR
        r_RR[RR] = r_s[isRR_s].mean()
    count_RR = (RR_s[:,None]==range(specrows)).sum(axis=0)

    if maskid == 'P000000N99':              
        RRaxis = np.argmin((raxis-r_RR)**2)
    elif maskid == 'P000000N03': 
        RRaxis = np.argmax(count_RR)

#   cull weak lines
    ptile = 100.*min(1.,5.*maskCols/count_RR.max())  # want like 5 brightest lines
    for RR in range(specrows):
        isRR_s = ok_s & np.in1d(bin_s,np.arange(bin0_RR[RR],bin1_RR[RR]+1))
        fluxmin = np.percentile(sex_js[2,isRR_s],100.-ptile)
        ok_s[isRR_s] &= (sex_js[2,isRR_s] > fluxmin) 
     
#   identify with mask rows R (assuming no gaps)
    RR_m = R_m + RRaxis - Raxis

#   find approximate grating shift in dispersion direction by looking for most common id error
    histc_b = np.zeros(60)
    for RR in range(specrows):
        isRR_s = ((RR_s==RR) & ok_s)
        cerr_MS = (c_s[None,isRR_s] - col_ml[RR_m==RR].ravel()[:,None]) 
        histc_b += np.histogram(cerr_MS.ravel(),bins=60,range=(-150,150))[0]
    cshift = 5*np.argmax(histc_b) - 150
    col_ml += cshift   

#   identify wavelength and mask column with spots in each spectrum
    isfound_s = np.zeros((spots),dtype=bool)
    bintol = 16/cbin                                  # 2 arcsec tolerance for line ID
    R_s = -np.ones(spots,dtype=int)
    C_s = -np.ones(spots,dtype=int)
    l_s = -np.ones(spots,dtype=int)
    m_s = -np.ones(spots,dtype=int)
    cerr_s = np.zeros(spots)
    rmscol = 0.

    for RR in range(specrows):      # _S spot in spectrum, _P (mask column, line)
        isRR_m = (RR_m==RR) 
        isRR_s = ((RR_s==RR) & ok_s)
        cerr_PS = (c_s[None,isRR_s] - col_ml[isRR_m].ravel()[:,None])
        Spots = isRR_s.sum()
        Possibles = col_ml[isRR_m].size
        Cols = Possibles/lines
        P_S = np.argmin(np.abs(cerr_PS),axis=0)
        cerr_S = cerr_PS[P_S,range(isRR_s.sum())]
        isfound_S = (np.abs(cerr_S) < bintol)
        M_P,l_P = np.unravel_index(np.arange(Possibles),(Cols,lines))
        m_P = np.where(isRR_m)[0][M_P]
        m_S = m_P[P_S]
        C_P = C_m[m_P]
        C_S = C_P[P_S]
        l_S = l_P[P_S]
        s_S = np.where(isRR_s)[0]
        R_s[isRR_s] = RR + Raxis-RRaxis
        cerr_s[s_S] = cerr_S
        C_s[s_S[isfound_S]] = C_S[isfound_S]
        l_s[s_S[isfound_S]] = l_S[isfound_S]
        m_s[s_S[isfound_S]] = m_S[isfound_S]
        isfound_s[s_S] |= isfound_S
        rmscol += (cerr_S[isfound_S]**2).sum()

#   cull wavelengths to _L with < 1/2 Mask Rows or Cols 
    ok_s &= isfound_s
    ok_l = np.zeros((lines),dtype=bool)
    for line in range(lines):
        lRows = np.unique(R_s[l_s==line]).shape[0]
        lCols = np.unique(C_s[l_s==line]).shape[0]
        ok_l[line] = ((lRows>=maskRows/2) & (lCols>=maskCols/2))
    l_L = np.where(ok_l)[0]
    wav_L = wav_l[l_L]
    Lines = l_L.shape[0]

    ok_s &= np.in1d(l_s,l_L)

#   tabulate good catalog spots (final _S)
    s_S = np.where(ok_s)[0]
    r_S = r_s[s_S]
    c_S = c_s[s_S]
    cerr_S = cerr_s[s_S]
    R_S = R_s[s_S]
    C_S = C_s[s_S]
    l_S = l_s[s_S]
    Spots = ok_s.sum()

    rshift = r_S[R_S==Raxis].mean() - raxis
    cshift += (c_S - col_ml[m_s[s_S],l_S]).mean()
    rmscol = np.sqrt(rmscol/Spots)

    np.savetxt("cat_S.txt",np.vstack((s_S,r_S,c_S,R_S,C_S,l_S,cerr_S)).T,  \
        fmt="%5i %8.2f %8.2f %5i %5i %5i %8.2f") 

    print "\nSpec_Spots Lines rshift  cshift      rms\n                pixels   pixels   pixels"
    print "  %5i %5i %8.1f %8.1f %8.1f" % (Spots,np.unique(l_S).shape[0],rshift,cshift,rmscol)
    print "\nLineno   Wavel   spots  Rows  Cols"
    for L in range(Lines):
        line = l_L[L]
        lRows = np.unique(R_S[l_S==line]).shape[0]
        lCols = np.unique(C_S[l_S==line]).shape[0]
        lspots = (l_S==line).sum()
        print " %5i %8.2f %5i %5i %5i" % (line,wav_l[line],lspots,lRows,lCols)

    sexcols = sex_js.shape[0]
    sexdata_jfS = np.zeros((sexcols,flexposns,Spots))
    sexdata_jfS[:,catpos] = sex_js[:,ok_s]
    xcenter_L = col_ml[maxis,l_L]
    ycenter = raxis + rshift 
    if option == "filesave":
        np.savetxt(prefix+"Spots.txt",sexdata_jfS[:,catpos].T,   \
            fmt=2*"%9.2f "+"%9.0f "+"%9.1f "+"%4i "+"%6.2f "+3*"%7.2f "+"%11.3e")

#   find spots in flexure series, in order of increasing abs(rho), and store sextractor output

    row_fLd = np.zeros((flexposns,Lines,2))
    col_fLd = np.zeros((flexposns,Lines,2))

    print "\n     fits      rho  line  spots  rshift  cshift  rslope  cslope  rmserr "
    print   "               deg   Ang         arcsec  arcsec  arcmin  arcmin   bins"

    for dirn in (1,-1):
        refpos = catpos
        posdirlist = np.argsort(dirn*rho_f)
        poslist = posdirlist[dirn*rho_f[posdirlist] > rho_f[refpos]]

        for fpos in poslist:
            col_S,row_S = sexdata_jfS[0:2,refpos,:]   
            sex_js = sextract(fitslist[fpos],"sexwt.fits",deblend=deblend)     

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

#            np.savetxt("motion_"+str(fpos)+".txt",np.vstack((isfound_S,l_S,drow_S,dcol_S)).T,fmt="%3i %3i %8.2f %8.2f")

        # Compute flexure image motion parameters for each line
            for L in range(Lines):
                ok_S = ((l_S == l_L[L]) & isfound_S)
                row_fLd[fpos,L],rowchi,d,d,d = \
                  np.polyfit(sexdata_jfS[0,catpos,ok_S]-xcenter_L[L],drow_S[ok_S],deg=1,full=True)
                col_fLd[fpos,L],colchi,d,d,d = \
                  np.polyfit(sexdata_jfS[1,catpos,ok_S]-ycenter,dcol_S[ok_S],deg=1,full=True)
                rms = np.sqrt((rowchi+colchi)/(2*ok_S.sum()))

                print ("%12s %5.0f %5i %5i "+5*"%7.2f ") % (image_f[fpos], rho_f[fpos], wav_L[L],  \
                    ok_S.sum(),row_fLd[fpos,L,1]*rbin*pix_scale, col_fLd[fpos,L,1]*cbin*pix_scale, \
                    60.*np.degrees(row_fLd[fpos,L,0]),-60.*np.degrees(col_fLd[fpos,L,0]), rms)
                if option == "filesave":
                    np.savetxt(prefix+"flex_"+str(fpos)+".txt",np.vstack((isfound_S,drow_S,dcol_S)).T,  \
                        fmt = "%2i %8.3f %8.3f")
                    np.savetxt(prefix+"sextr_"+str(fpos)+".txt",sexdata_jfS[:,fpos].T)
            print 

#   make plots

    fig,plot_s = plt.subplots(2,1,sharex=True)

    plt.xlabel('Rho (deg)')
    plt.xlim(-120,120)
    plt.xticks(range(-120,120,30)) 
    fig.set_size_inches((8.5,11))
    fig.subplots_adjust(left=0.175)

    plot_s[0].set_title(str(dateobs)+[" Imaging"," Spectral"][isspec]+" Flexure") 
    plot_s[0].set_ylabel('Mean Position (arcsec)')
    plot_s[0].set_ylim(-0.5,4.)
    plot_s[1].set_ylabel('Rotation (arcmin ccw)')      
    plot_s[1].set_ylim(-10.,6.)

    lbl_L = [("%5.0f") % (wav_L[L]) for L in range(Lines)]
    color_L = 'bgrcmykw'
    for L in range(Lines):
      plot_s[0].plot(rho_f,row_fLd[:,L,1]*rbin*pix_scale,   \
            color=color_L[L],marker='D',markersize=8,label='row '+lbl_L[L])
      plot_s[1].plot(rho_f,60.*np.degrees(row_fLd[:,L,0]),  
            color=color_L[L],marker='D',markersize=8,label='row '+lbl_L[L])
    collbl = 'col'+lbl_L[0]
    for L in range(Lines): 
      plot_s[0].plot(rho_f,col_fLd[:,L,1]*cbin*pix_scale,   \
            color=color_L[L],marker='s',markersize=8,label=collbl)
      plot_s[1].plot(rho_f,-60.*np.degrees(col_fLd[:,L,0]), \
            color=color_L[L],marker='s',markersize=8,label=collbl)
      collbl = ''
    plot_s[0].legend(fontsize='medium',loc='upper center')                     
    plotfile = str(dateobs)+['_imflex.pdf','_grflex.pdf'][isspec]
    plt.savefig(plotfile,orientation='portrait')

    if os.name=='posix':
        if os.popen('ps -C evince -f').read().count(plotfile)==0: os.system('evince '+plotfile+' &')

    os.remove("out.txt")
    os.remove("qred_thrufoc.param") 
    os.remove("sexwt.fits")                           
    return

# ---------------------------------------------------------------------------------

if __name__=='__main__':
    imagefits=sys.argv[1]
    fitslist=sys.argv[2:]
    option = ""
    if (fitslist[-1] == "filesave"):
        option = fitslist.pop()
    flexure_rssspec(imagefits,fitslist,option)

# debug:
# cd /d/pfis/khn/20110804/sci
# flexure_rssspec.py m*E304.fits m*E4??.fits
# cd /d/pfis/khn/20160725/sci
# flexure_rssspec.py m*0063.fits m*006[8-9].fits m*007[0-6].fits filesave
