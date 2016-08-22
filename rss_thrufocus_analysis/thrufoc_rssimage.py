#! /usr/bin/env python

# Compute RSS focus plane for thrufocus run, any imaging mask

import os, sys, time
import numpy as np
from datetime import datetime
from scipy import interpolate as ip
from scipy import linalg as la
import pyfits

np.set_printoptions(threshold=np.nan)

datadir = "/d/carol/Synched/software/SALT/pipeline/"    

# ---------------------------------------------------------------------------------
def sextract(fits,binmask="",deblend=.02):
    global gap_gc, rfovpix
# run sextractor to find objects.  Altered parameter list qred_thrufoc.param in cwd
    global hdulist_image, shape_d
    hdulist_image = pyfits.open(fits)
    cbin,rbin = map(int,(hdulist_image[0].header["CCDSUM"].split(" ")))
    maskid = hdulist_image[0].header["MASKID"]
    isspec = hdulist_image[0].header["GR-STATE"][1] == "4"
    image_rc = np.copy(hdulist_image["SCI"].data)
    rows,cols = image_rc.shape
    sigma = 6.                    
    pix_scale=0.125 
    sat=0.99*image_rc.max()
    if isspec: rfovpix = 4000.              # fov larger than detector
    else: rfovpix = 250./pix_scale

# if no binmask specified, mask out gaps and area outside fov.  Coordinates are bins relative to ccd center, assumed to be center of image
    if binmask == "":
        binmask = "sexwt.fits"
        left,right = np.array([np.where(image_rc>0)[1].min(), np.where(image_rc>0)[1].max()])
        gap_gc = np.array([[left+2048/cbin-2,cols/2-1024/cbin+2],[cols/2+1024/cbin-2,right-2048/cbin+2]]) - cols/2
        row_r = np.arange(rows) - rows/2
        col_c = np.arange(cols) - cols/2
        mask_rc = (np.logical_not(((col_c>gap_gc[0,0]) & (col_c<gap_gc[0,1]))[None,:] | \
                              ((col_c>gap_gc[1,0]) & (col_c<gap_gc[1,1]))[None,:] | \
                              (np.sqrt((rbin*row_r[:,None])**2+(cbin*col_c[None,:])**2) > rfovpix))).astype(np.uint8)
        pyfits.PrimaryHDU(mask_rc).writeto(binmask, clobber = "True")
 
    cmd= ('sex %s -c '+datadir+'qred.sex -PARAMETERS_NAME %s -WEIGHT_IMAGE %s -WEIGHT_TYPE MAP_WEIGHT -CATALOG_NAME %s \
        -DETECT_THRESH %f -SATUR_LEVEL %f -DEBLEND_MINCONT %f') % (fits, "qred_thrufoc.param", binmask, "out.txt", sigma, sat, deblend)
    os.system(cmd+" &> /dev/null")
#    os.system(cmd)
    sex_js = np.loadtxt("out.txt").transpose()
    return sex_js

# ---------------------------------------------------------------------------------
def thrufoc_rssimage(fitslist,option=""):
    global sex_js, rd, B, pixarcsec, gridsize, niter
    global focus_f, rc0_dgg, usespots_gg, fwhminterp_g                                            # for cube analysis

    if option == "filesave": prefix = raw_input("file prefix: ")
    pixel = 15.                                                                  # pixel size in microns
    pix_scale=0.125 
    turnfac_d = [8.0,6.7]                                                        # arcmin/turn for tip/tilt adustment
    sexparams = ["X_IMAGE","Y_IMAGE","FLUX_ISO","FLUX_MAX","FLAGS","CLASS_STAR",    \
                "X2WIN_IMAGE","Y2WIN_IMAGE","XYWIN_IMAGE","ERRX2WIN_IMAGE"]
    np.savetxt("qred_thrufoc.param",sexparams,fmt="%s")
    fmaxcol,flagcol,xvarcol,yvarcol,xerrcol = (3,4,6,7,9)                        # column nos (from 0) of data in sextractor
    imagestooclosefactor = 3.0                                                   # too close if factor*sep < sqrt(var)
    gaptooclose = 1.25                                                           # arcsec
    edgetooclose = 1.25                                                          # arcsec
    rattolerance = 0.25    
    toofaint = 250.                                                              # FMAX counts
    galaxydelta = 0.4                                                            # arcsec
    MOSimagelimit = 1.                                                           # arcsec

    focposns = len(fitslist)
    catposn = focposns/2
    hdr = pyfits.getheader(fitslist[catposn],0)
    cbin,rbin = map(int,(hdr["CCDSUM"].split(" ")))
    rows,cols = pyfits.getdata(fitslist[catposn]).shape
    maskid =  (hdr["MASKID"]).strip() 
    filter =  (hdr["FILTER"]).strip()
    if "COLTEM" in hdr:                                                         # allow for annoying SALT fits version change
        coltem = float(hdr["COLTEM"]) 
        camtem = float(hdr["CAMTEM"])
    else:
        coltem = float(hdr["COLTEMP"]) 
        camtem = float(hdr["CAMTEMP"])     
    focuscat = hdr["FO-POS-V"]

    print str(datetime.now()), "\n"
    if ((maskid[0:2]=="P0") & (maskid[-2:]!="00")):
        deblend = 1      # MOS or pinhole mask - no deblending
        print "Mask:      ", maskid," - MOS or pinhole mask - no deblending"
    else:
        deblend=.02      # higher than default .005 to prevent deblending of out-of-focus stars
        print "Mask:      ", maskid," - Imaging mask - minimal deblending"
    print "Filter:    ", filter
    print "Col Temp:   %7.1f" % (coltem)
    print "Cam Temp:   %7.1f" % (camtem)

#   make catalog of stars using middle image (capital S)
    sex_jS = sextract(fitslist[catposn],deblend=deblend)
    sexcols,stars = sex_jS.shape    
    col_S,row_S = sex_jS[0:2,:]

#    for S in range(stars): print (2*"%9.2f "+"%9.0f "+"%9.1f "+"%4i "+"%6.2f "+3*"%7.2f "+"%11.3e ") % tuple(sex_jS[:,S])

#   find stars in focus series, and store sextractor output
    image_f = []
    focus_f = np.zeros((focposns))
    sexdata_jfS = np.zeros((sexcols,focposns,stars))
    print "\n fits  focus  stars  rshift  cshift   xfwhm   yfwhm   fwhm   rmserr "
    print   "       (micr)        arcsec  arcsec  arcsec  arcsec  arcsec   bins"
    for fpos in range(focposns):   
        sex_js = sextract(fitslist[fpos],"sexwt.fits",deblend=deblend)     # use standard mask made for middle image
        image_f.append(fitslist[fpos].split(".")[0][-6:])
        bintol = 32/cbin                                                   # 4 arcsec tolerance for finding star
        binsqerr_sS = (sex_js[1,:,None] - row_S[None,:])**2 + (sex_js[0,:,None] - col_S[None,:])**2
        S_s = np.argmin(binsqerr_sS,axis=1)
    # First compute image shift by averaging small errors
        rowerr_s = sex_js[1] - row_S[S_s]
        colerr_s = sex_js[0] - col_S[S_s]
        hist_r,bin_r = np.histogram(rowerr_s,bins=32,range=(-2*bintol,2*bintol))    
        drow = rowerr_s[(rowerr_s > bin_r[np.argmax(hist_r)]-bintol) & (rowerr_s < bin_r[np.argmax(hist_r)]+bintol)].mean()
        hist_c,bin_c = np.histogram(colerr_s,bins=32,range=(-2*bintol,2*bintol))    
        dcol = colerr_s[(colerr_s > bin_r[np.argmax(hist_r)]-bintol) & (colerr_s < bin_r[np.argmax(hist_r)]+bintol)].mean()
    # Now find the closest ID
        binsqerr_sS = (sex_js[1,:,None] - row_S[None,:] -drow)**2 + (sex_js[0,:,None] - col_S[None,:] -dcol)**2
        binsqerr_s = binsqerr_sS.min(axis=1)
        isfound_s = binsqerr_s < bintol**2
        S_s = np.argmin(binsqerr_sS,axis=1)
        isfound_s &= (binsqerr_s == binsqerr_sS[:,S_s].min(axis=0))
        sumsq = binsqerr_sS[isfound_s,S_s[isfound_s]].sum()
        rms = np.sqrt(sumsq/isfound_s.shape[0])
        focus_f[fpos] = pyfits.getheader(fitslist[fpos],0)["FO-POS-V"]
        sexdata_jfS[:,fpos,S_s[isfound_s]] = sex_js[:,isfound_s]
        xwidth = (cbin*pix_scale)*2.*np.median(np.sqrt(sex_js[xvarcol,isfound_s]))                # in arcsec
        ywidth = (rbin*pix_scale)*2.*np.median(np.sqrt(sex_js[yvarcol,isfound_s]))
        fwhm = np.sqrt(xwidth**2+ywidth**2)
        print ("%6s %5.0f %5i "+6*"%7.2f ") % \
            (image_f[fpos], focus_f[fpos], isfound_s.sum(), drow*rbin*pix_scale, dcol*cbin*pix_scale, xwidth, ywidth, fwhm, rms)

#   assemble good focus curves, with enough focus samples of similar-looking images

    if option == "filesave": 
        np.savetxt(prefix+'_fmax.txt',sexdata_jfS[fmaxcol].T,fmt="%9.1f")
        np.savetxt(prefix+'_xvar.txt',sexdata_jfS[xvarcol].T,fmt="%7.2f")
        np.savetxt(prefix+'_yvar.txt',sexdata_jfS[yvarcol].T,fmt="%7.2f")

    use_fS = sexdata_jfS[0] > 0.
    cullflag_fS = np.logical_not(use_fS).astype(int)                                # cullflag=1 for object not found

#   throw out images that are too close
    xpos_fS,ypos_fS = sexdata_jfS[0:2]
    xvar_fS,yvar_fS = sexdata_jfS[xvarcol:xvarcol+2]
    argminsep_fS = np.argsort((xpos_fS[:,:,None]-xpos_fS[:,None,:])**2 + (ypos_fS[:,:,None]-ypos_fS[:,None,:])**2)[:,:,1]
    nbrdata_jfS = np.zeros_like(sexdata_jfS)
    for f in range(focposns): nbrdata_jfS[:,f,:] = sexdata_jfS[:,f,argminsep_fS[f]] 
    nbrxpos_fS,nbrypos_fS = nbrdata_jfS[0:2]
    nbrxvar_fS,nbryvar_fS = nbrdata_jfS[xvarcol:xvarcol+2]
    dxmin_fS = nbrxpos_fS - xpos_fS; dymin_fS = nbrypos_fS - ypos_fS
    drmin_fS = np.zeros_like(xpos_fS); minseprat_fS = np.zeros_like(xpos_fS)
    drmin_fS[use_fS] = np.sqrt(dxmin_fS[use_fS]**2 + dymin_fS[use_fS]**2)
    minseprat_fS[use_fS] = drmin_fS[use_fS] \
       / (np.sqrt(xvar_fS[use_fS]*(dxmin_fS[use_fS]/drmin_fS[use_fS])**2 + yvar_fS[use_fS]*(dymin_fS[use_fS]/drmin_fS[use_fS])**2)  \
       + np.sqrt(nbrxvar_fS[use_fS]*(dxmin_fS[use_fS]/drmin_fS[use_fS])**2 + nbryvar_fS[use_fS]*(dymin_fS[use_fS]/drmin_fS[use_fS])**2))

    print "\nImages remaining after culling for:          Final focus curves: \n begin sep  flux shape  edge count blobs"
    print "%4i " % use_fS.sum(),

    ok_fS = (minseprat_fS > imagestooclosefactor)
    cullflag_fS += 2*(use_fS & np.logical_not(ok_fS)).astype(int)                       # cullflag = 2 for too close
    use_fS &= ok_fS

    print "%4i " % use_fS.sum(),

#   throw out images that are too faint
    ok_fS = (sexdata_jfS[fmaxcol] > toofaint)                               
    cullflag_fS += 3*(use_fS & np.logical_not(ok_fS)).astype(int)                       # cullflag = 3 for too faint
    use_fS &= ok_fS
    print "%4i " % use_fS.sum(),                                                                          

#   For images much larger than spectrograph psf, throw out images with unusual shapes.  
#       Then use xvar instead of total FWHM if yvar/xvar > 4 (MOS mask).
    fwhm_fS = 2.*np.sqrt((cbin*pix_scale)**2*xvar_fS + (rbin*pix_scale)**2*yvar_fS)
    if np.median(fwhm_fS[use_fS]) > MOSimagelimit :
        YXrat_fS = np.zeros_like(xvar_fS)
        YXrat_fS[use_fS] = (yvar_fS[use_fS]/xvar_fS[use_fS])*(rbin/cbin)**2
        YXrat = np.median(YXrat_fS[use_fS])                                                 
        ok_fS = (YXrat_fS > (1.-rattolerance)*YXrat) & (YXrat_fS < (1.+rattolerance)*YXrat)
        cullflag_fS += 4*(use_fS & np.logical_not(ok_fS)).astype(int)                    # cullflag = 4 for bad shape
        use_fS &= ok_fS   
        YXrat = np.median(YXrat_fS[use_fS])                                              # median shape of good images
        if YXrat > 4.:
            fwhm_fS = 2.*np.sqrt((cbin*pix_scale)**2*xvar_fS)
            print "\nDetect MOS mask: Y/X variance ratio = %8.2f" % (YXrat)
 
    print "%4i " % use_fS.sum(),

#   throw out images too near gaps and edge
    ok_fS = (np.absolute(xpos_fS[:,:,None] - (gap_gc.flatten()[None,None,:] + cols/2)).min(axis=2) > gaptooclose)
    ok_fS &= (np.sqrt((rbin*(ypos_fS-rows/2.))**2+(cbin*(xpos_fS-cols/2.))**2) < (rfovpix - edgetooclose/pix_scale))
    cullflag_fS += 5*(use_fS & np.logical_not(ok_fS)).astype(int)                        # cullflag = 5 for too near edges
    use_fS &= ok_fS
    print "%4i " % use_fS.sum(),

#   throw out stars missing too many good images
    focuscount_S = use_fS.sum(axis=0)
    ok_S = (focuscount_S > focposns/2)
    cullflag_fS += 6*(use_fS & np.logical_not(ok_S[None,:])).astype(int)                 # cullflag = 6 for too few points in focus curve
    use_fS &= ok_S[None,:]
    use_S = use_fS.sum(axis=0) > 0
    print "%4i " % use_fS.sum(),         
                                                                                    
#   throw out "stars" with consistently large size (eg, galaxies)
    dfwhm_fS = np.zeros((focposns,stars));  offset_S = np.zeros(stars)
    for f in range(focposns): dfwhm_fS[f,use_fS[f]] = fwhm_fS[f,use_fS[f]] - np.median(fwhm_fS[f,use_fS[f]])
    for S in np.where(use_S)[0]: offset_S[S] = np.min(dfwhm_fS[use_fS[:,S],S])
    ok_S = (offset_S < galaxydelta)
    cullflag_fS += 7*(use_fS & np.logical_not(ok_S[None,:])).astype(int)                  # cullflag = 7 for blobs
    use_fS &= ok_S[None,:]
    use_f = use_fS.sum(axis=1) > 0
    use_S = use_fS.sum(axis=0) > 0
    print "%4i " % use_fS.sum(), "     %4i" % use_S.sum(), "\n"

    if option=="filesave":
        cullfile = open(prefix+"_cull.txt","w")
        print >> cullfile, "star   row   column   cullflag"
        for S in range(stars): 
            print  >> cullfile, ("%3i %7.1f %7.1f " % (S,row_S[S],col_S[S])), (focposns*"%2i " % tuple(cullflag_fS[:,S].astype(int)))

#   make interpolated master focus curve and find best focus shift for each star

    focus0 = focus_f[np.where(use_f)[0][0]]
    focus1 = focus_f[np.where(use_f)[0][-1]]
    dfocus = int(focus1 - focus0)
    focusgrid_F = np.arange(focus0,focus1+1,1.)
    fwhm_f = np.zeros(focposns)
    for f in range(focposns): fwhm_f[f] = np.median(fwhm_fS[f,use_fS[f]])
    minfocus = focus_f[np.argmin(fwhm_f)]
    focusrange = (np.argmin(fwhm_f)==0) | (np.argmin(fwhm_f)==fwhm_f.shape[0]-1)
    fwhm_F = ip.interp1d(focus_f,fwhm_f,kind ='cubic')(focusgrid_F)

    focshift_S = np.zeros(stars,dtype="int");    fwhmmin_S = np.zeros(stars) 

    for S in np.nonzero(use_S)[0]:
        focusS_f = focus_f[use_fS[:,S]]
        dfocusS = int(focusS_f[-1] - focusS_f[0])
        dfocusS0 = int(focusS_f[0] - focus0)
        minguess_f = fwhm_fS[np.where(use_fS[:,S])[0][1:-1],S]          # ignore 1st and last focal positions of this curve
        minfocusS = focusS_f[np.argmin(minguess_f)+1]                   # 1st guess for best focus to center search
        if (focusrange): minfocusS = minfocus                           # if median curve has no min, use default search
        focusgridS_F = np.arange(focusS_f[0],focusS_f[-1]+1,1.)
        fwhmS_F = ip.interp1d(focus_f[use_fS[:,S]],fwhm_fS[use_fS[:,S],S],kind ='cubic')(focusgridS_F)

        focshift_t = np.arange(max(-0.5*dfocus,minfocusS-minfocus-150), min(0.5*dfocus,minfocusS-minfocus+150)).astype(int)
        rms_t = np.zeros(focshift_t.shape[0])
        for t in range(focshift_t.shape[0]):
            idxS0 = np.minimum(dfocusS,np.maximum(0,focshift_t[t]-dfocusS0))
            idxS1 = np.minimum(dfocusS,np.maximum(0,dfocusS+focshift_t[t]-dfocusS0)) + 1
            idx0 = idxS0 - focshift_t[t] + dfocusS0
            idx1 = idxS1 - focshift_t[t] + dfocusS0 
            shaperr = fwhmS_F[idxS0:idxS1].shape[0] -  fwhm_F[idx0:idx1].shape[0]
            rms_t[t] = np.sqrt(((fwhmS_F[idxS0:idxS1] - fwhm_F[idx0:idx1])**2).mean())

        focshift_S[S] = focshift_t[np.argmin(rms_t)]
        fwhmmin_S[S] = fwhmS_F.min()   
        if S==4:
            print minfocus
            print focusS_f
            print dfocus
            print dfocusS0
            print minguess_f
            print minfocusS
            print focusgridS_F
            print focshift_t
            print rms_t
            print np.argmin(rms_t)
            exit()

    if option=="filesave":
        focplanefile = open(prefix+"_focplane.txt","w")
        print >> focplanefile, " star  col    row  focshift fwhmmin"
        print >> focplanefile, "       bin    bin    micron  arcsec"
        for S in np.where(use_S)[0]: print >> focplanefile, \
            ("%3i %7.1f %7.1f %7i %6.3f" % (S, col_S[S], row_S[S], focshift_S[S], fwhmmin_S[S]))

#   fit 1st order 2D polynomial to focshift (row,col)
    b_s = focshift_S[use_S]
    a_sd = np.vstack((np.ones_like(b_s),row_S[use_S]-rows/2,col_S[use_S]-cols/2)).T
    cof_d,sumsqerr = la.lstsq(a_sd,b_s)[0:2]
    sigma = np.sqrt(sumsqerr/b_s.shape[0])    
    alpha_dd = (a_sd[:,:,None]*a_sd[:,None,:]).sum(axis=0)
    eps_dd = la.inv(alpha_dd)                                        # errors from Bevington, p121
    err_d = sigma*np.sqrt(np.diagonal(eps_dd))
    tip,tilt = 60.*np.degrees(cof_d[1:3])/(pixel*rbin)
    tiperr,tilterr = 60.*np.degrees(err_d[1:3])/(pixel*rbin)

    print "Focal plane fit, err: ",
    if focusrange: print "   (Warning: mean focus is at limit)",
    print "\n       mean     tip     tilt        rms  \n      micron  arcmin  arcmin      micron"
    print "     %7.1f %7.2f %7.2f %10.1f"  % (cof_d[0]+minfocus,tip,tilt,sigma)
    print " +/- %7.1f %7.2f %7.2f"  % (err_d[0],tiperr,tilterr)

#   fit 2nd order 2D polynomial to focshift (row,col)
    aa_sd = np.vstack(((row_S[use_S]-rows/2)**2,(col_S[use_S]-cols/2)**2,(row_S[use_S]-rows/2)*(col_S[use_S]-cols/2))).T
    a_sd = np.hstack((a_sd,aa_sd))
    cof_d,sumsqerr = la.lstsq(a_sd,b_s)[0:2]
    sigma = np.sqrt(sumsqerr/b_s.shape[0])
    alpha_dd = (a_sd[:,:,None]*a_sd[:,None,:]).sum(axis=0)
    eps_dd = la.inv(alpha_dd)                                        
    err_d = sigma*np.sqrt(np.diagonal(eps_dd))    
    curv = 0.5*(cof_d[3]*(rfovpix/rbin)**2 + cof_d[4]*(rfovpix/cbin)**2)
    curverr = 0.5*(err_d[3]*(rfovpix/rbin)**2 + err_d[4]*(rfovpix/cbin)**2)

    print "\nFit focal surface to quadratic:"
    print " curvature: ", "%7.1f +/- %7.1f" % (curv,curverr)

    phirim_p = np.arange(360.)
    x_p = (rfovpix/cbin)*np.sin(np.radians(phirim_p))
    y_p = (rfovpix/rbin)*np.cos(np.radians(phirim_p))
    b_s = cof_d[0] + cof_d[1]*y_p + cof_d[2]*x_p + cof_d[3]*y_p**2 + cof_d[4]*x_p**2  + cof_d[5]*x_p*y_p
    a_sd = np.vstack((np.ones_like(b_s),y_p,x_p)).T
    cof_d = la.lstsq(a_sd,b_s)[0]
    tip,tilt = 60.*np.degrees(cof_d[1:3])/(pixel*rbin)

    print "Rim focal plane fit, err: ",
    print "\n       mean     tip     tilt  \n      micron  arcmin  arcmin"
    print "     %7.1f %7.2f %7.2f"  % (cof_d[0]+minfocus,tip,tilt)

#   make output fits cube: tile(64x64 star image vs focus) vs star
#   center each tile (to closest pixel) on star position at closest to best focus
    if option=="filesave":
        cubeout_Sfcr = np.zeros((stars,focposns,64/cbin,64/rbin),dtype="float32")
        fidx_S = -np.ones((stars))
        fmin_S = minfocus + focshift_S
        hduoutl = pyfits.open(fitslist[catposn])
        for f in range(focposns):   
            cubeout_cr = pyfits.open(fitslist[f])["SCI"].data.T
            for S in np.where(use_S)[0]:
                fidx_S[S] = np.argmin(np.abs(fmin_S[S]-focus_f))
                col_S[S] = sexdata_jfS[0,fidx_S[S],S].round()
                row_S[S] = sexdata_jfS[1,fidx_S[S],S].round()
                cubeout_Sfcr[S,f] = \
                  cubeout_cr[col_S[S]-32/cbin:col_S[S]+32/cbin,row_S[S]-32/rbin:row_S[S]+32/rbin]
        print minfocus
        print focshift_S
        print fmin_S
        hduoutl["SCI"].data = cubeout_Sfcr.reshape((stars,-1,64/rbin)).transpose((0,2,1))
        hduoutl[0].header.update("FOCUS",(focposns*"%i " % tuple(focus_f.astype(int))),'')
        c1 = pyfits.Column(name='bestfoc',format='J',array=fmin_S)
        fminhdu = pyfits.new_table([c1,])
        fminhdu.update_ext_name("BESTFOC")
        hduoutl.append(fminhdu)  
        hduoutl.writeto(prefix+"_cube.fits", clobber = "True")

        cubedatafile = open(prefix+"_cubedata.txt","w")
        print >> cubedatafile, "  star  col   row   focus  img"
        print >> cubedatafile, "        bin   bin  micron"
        for S in np.where(use_S)[0]: print >> cubedatafile, \
            ("%3i %6.0f %6.0f %6.0f %3i" % (S, col_S[S], row_S[S], minfocus+focshift_S[S], fidx_S[S]))

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

    thrufoc_rssimage(fitslist,option)
