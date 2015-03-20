#! /usr/bin/env python

# Compute RSS imaging flatfield for a specific image

import os, sys, time
from datetime import datetime
import numpy as np
from scipy.ndimage import interpolation as nd
from scipy import interpolate as ip
import pyfits
from math import pi

np.set_printoptions(threshold=np.nan)

trkgrid = "1503"                                        # yymm version of flat-track grid                       
datadir = "/d/carol/Synched/software/SALT/pipeline"
distfile = datadir+"/spectrograph/imgdist.txt"
trackfile = datadir+"/telescope/track_"+trkgrid+"_telflat.fits"

# ---------------------------------------------------------------------------------
def img_recenter(ar,bincen_d):
# move bincen to center of array
    shape_d = np.array(ar.shape)
    rows,cols = shape_d
    roff,coff = (shape_d/2 - bincen_d).astype(int)
    arr = np.zeros_like(ar)
    arr[max(0,roff):min(rows,rows+roff),max(0,coff):min(cols,cols+coff)] \
          = ar[max(0,-roff):min(rows,rows-roff), max(0,-coff):min(cols,cols-coff)] 
    return arr

# ---------------------------------------------------------------------------------
def blksmooth2d(ar,blk,maskcol_c):
# blkaverage (ignoring mask columns, plus blk values < 1/3 max), then spline interpolate result
    okbin = ar>0.
    gridyx = np.indices(ar.shape)
    blks_d = (np.asarray(ar.shape)/blk).astype(int)

    minlim = 0.33
    arr = ar; arr[:,maskcol_c] = 0.
    smbin = okbin; smbin[:,maskcol_c] = 0.
    offset_d = (np.asarray(ar.shape) - blks_d*blk)/2
    arrr = arr[offset_d[0]:offset_d[0]+blks_d[0]*blk, offset_d[1]:offset_d[1]+blks_d[1]*blk]    # cut to integral blocks
    sum = arrr.reshape(blks_d[0], blk, blks_d[1], blk).sum(axis=3).sum(axis=1)
    count = (arrr.reshape(blks_d[0], blk, blks_d[1], blk) > 0).sum(axis=3).sum(axis=1)
    arblk = np.zeros_like(sum)
    okblk = count > (blk**2)/4
    arblk[okblk] = sum[okblk]/count[okblk]
    okblk &= arblk > minlim*arblk.max()
    griddy,griddx = (gridyx*smbin.astype(int))[:,offset_d[0]:offset_d[0]+blks_d[0]*blk, offset_d[1]:offset_d[1]+blks_d[1]*blk]
    gridblky,gridblkx = np.zeros_like(sum),np.zeros_like(sum)
    gridblky[okblk] = (griddy.reshape(blks_d[0], blk, blks_d[1], blk).sum(axis=3).sum(axis=1))[okblk]/count[okblk]
    gridblkx[okblk] = (griddx.reshape(blks_d[0], blk, blks_d[1], blk).sum(axis=3).sum(axis=1))[okblk]/count[okblk]
    arsmooth = np.nan_to_num(ip.griddata((gridblky[okblk],gridblkx[okblk]), arblk[okblk], (gridyx[0],gridyx[1]), method="linear"))

    return arsmooth

# ---------------------------------------------------------------------------------
def ampgains(image_rc,rcbin_d):
# find amplifier columns, avoiding amplifier boundaries for fit
    centercol_s = np.zeros(3,dtype=int)                                 # center column in each ccd (start of righthand amp)
    col1_a = np.zeros(6,dtype=int); col2_a = np.zeros(6,dtype=int); 
    ampcol_ac = [np.array(a) for a in ([0],[0],[0],[0],[0],[0])]        # columns in each amp (list of arrays)    
    gain_a = np.ones(6)
    shape_d = np.asarray(image_rc.shape)

    centercol_s[0] = np.where(image_rc>0)[1].min() + 1024/rcbin_d[1] - 1
    centercol_s[1] = shape_d[1]/2 
    centercol_s[2] = np.where(image_rc>0)[1].max() - 1024/rcbin_d[1] + 1

#   set columns to be used for gain determination (avoiding partial bins)
    ampwidth_a = 1024/rcbin_d[1] - [ 2,2,1,1,2,2 ]
    for amp in range(0,6,2):
        col1_a[amp] = centercol_s[amp/2] - ampwidth_a[amp]
        col2_a[amp] = centercol_s[amp/2] - 1
        col1_a[amp+1] = centercol_s[amp/2] + 1
        col2_a[amp+1] = centercol_s[amp/2] + ampwidth_a[amp]

# average together 256 pixels at center row, ignoring zeros
    rows = 256/rcbin_d[0]
    cols = 64/rcbin_d[1]
    xprof_rc = image_rc[(shape_d[0]-rows)/2:(shape_d[0]+rows)/2 , :]
    xprof_c = xprof_rc.sum(axis=0)
    count = (xprof_rc>0).sum(axis=0)
    okcol = np.where(count>0)
    xprof_c[okcol] /= count[okcol]

# find gain ratios by linear extrapolation of outer 64 pixels in each amplifier
    for amp in range(1,6):
        matchcol = (col2_a[amp-1] + col1_a[amp])/2. 
        left = np.polyval(np.polyfit(range(col2_a[amp-1]-cols,col2_a[amp-1]),xprof_c[col2_a[amp-1]-cols:col2_a[amp-1]],1),matchcol)
        right = np.polyval(np.polyfit(range(col1_a[amp],col1_a[amp]+cols),xprof_c[col1_a[amp]:col1_a[amp]+cols],1),matchcol)
        gain_a[amp] = gain_a[amp-1]*right/left
    gain_a /= gain_a[2]                 # define lefthand amplifier in central CCD as gain 1

# return column range for each gain, plus gap edges to avoid
    ampcol_ac[0] = np.arange(centercol_s[0])
    ampcol_ac[1] = np.arange(centercol_s[0],(centercol_s[0]+centercol_s[1])/2)  
    ampcol_ac[2] = np.arange((centercol_s[0]+centercol_s[1])/2,centercol_s[1]) 
    ampcol_ac[3] = np.arange(centercol_s[1],(centercol_s[1]+centercol_s[2])/2)  
    ampcol_ac[4] = np.arange((centercol_s[1]+centercol_s[2])/2,centercol_s[2])
    ampcol_ac[5] = np.arange(centercol_s[2], image_rc.shape[1])  
  
    gapedge1 = np.arange(col2_a[1]-16/rcbin_d[1], col2_a[1]+16/rcbin_d[1])   # gap1 area to avoid
    gapedge2 = np.arange(col1_a[4]-16/rcbin_d[1], col1_a[4]+16/rcbin_d[1])   # gap2 area to avoid
    maskcol_c = np.append(gapedge1,gapedge2)

    return gain_a,ampcol_ac,maskcol_c   

# ---------------------------------------------------------------------------------
def lampflat_rssimage(fits, gaincor=True):
    global gain_a,ampcol_ac,maskcol_c
    hdulist_image = pyfits.open(fits)
    flat_rc = np.copy(hdulist_image["SCI"].data)
    shape_d = np.asarray(flat_rc.shape)
    rcbin_d = np.array(hdulist_image[0].header["CCDSUM"].split(" ")).astype(int)
    dtype_image = flat_rc.dtype

# separate amplifiers and correct for gain.  If gaincor False, prefer gains from previous (sky) flat contribution
    if gaincor: gain_a,ampcol_ac,maskcol_c = ampgains(flat_rc,rcbin_d)
    for amp in range(6): flat_rc[:,ampcol_ac[amp]] /= gain_a[amp]

# remove low-frequency shape by blocking to 32 unbinned pixels, then fitting with spline
    blk = 32/rcbin_d.min()
    flatsmooth_rc = blksmooth2d(flat_rc,blk,maskcol_c)
    okpix = flatsmooth_rc > flatsmooth_rc.max()/4.
    flat_rc[okpix] /= flatsmooth_rc[okpix]
    badpix = flatsmooth_rc <= flatsmooth_rc.max()/4.
    flat_rc[badpix] = 0.

# put gain back in unless gaincor False
    if gaincor: 
        for amp in range(6): flat_rc[:,ampcol_ac[amp]] *= gain_a[amp] 
  
    return flat_rc

# ---------------------------------------------------------------------------------
def skyflat_rssimage(fits, gaincor=True):
#   gaincor:  if False, leave flat smooth (no gain correction)

    global gain_a,ampcol_ac,maskcol_c
    hdulist_image = pyfits.open(fits)
    flat_rc = np.copy(hdulist_image["SCI"].data)
    shape_d = np.asarray(flat_rc.shape)
    rcbin_d = np.array(hdulist_image[0].header["CCDSUM"].split(" ")).astype(int)
    dtype_image = flat_rc.dtype

# make sextractor weight file with flat < 60% max = 0
    maxval = flat_rc.max()
    hist_v, bin_v = np.histogram(flat_rc,bins=100,range=(0.01,0.99*maxval))  # get rid of zeroes, sat'n in histogram
    maxbin = np.where(hist_v > 1.e-4*hist_v.sum())[0].max()                 # get rid of cr's, bright stars in histogram
    hist_v, bin_v = np.histogram(flat_rc,bins=100,range=(0.01,bin_v[maxbin]))  
    maxhistarg = (hist_v*bin_v[:-1]).argmax() 
    maxhistval = hist_v[maxhistarg]
    maxflatval = bin_v[np.where(hist_v > 0.5*maxhistval)[0].max()]          # max flat value is where hist drops 2x from hist max

    flatwt_rc = (flat_rc > 0.6*maxflatval).astype(np.int16)
    hdulist_image["SCI"].data = flatwt_rc
    wtfile = fits.split("/")[-1].split(".")[0] + "_wt.fits"
    hdulist_image.writeto(wtfile, clobber = "True")                    

# run sextractor on it to find objects
    pix_scale=0.125
    r_ap=2.0/(pix_scale*rcbin_d.min())         
    sat=0.99*maxval
    Apsf,Bpsf = 1.280, 0.241                                                      # PSF wing calibration 11 Mar'15, with 20% boost for bright                              # PSF wing calibration
    psflim = 0.01                                                                 # limit star wings to this fraction of sky
    outtxt=(fits.split("/")[-1]).split('.fits')[0]+'.txt'
    cmd='sex %s -c $PYTHONPATH/qred.sex -CATALOG_NAME %s -WEIGHT_IMAGE %s -WEIGHT_TYPE MAP_WEIGHT -PHOT_APERTURES %f -SATUR_LEVEL %f &> \
        /dev/null' % (fits, outtxt, wtfile, r_ap, sat)
    os.system(cmd)

# remove objects from image.                                                              
    y_im,x_im = np.indices(shape_d)       
    x_s,y_s,flux_s,fluxmax_s,fwhm_s = np.loadtxt(outtxt,usecols=(0,1,4,7,10),unpack=True)
    fwhm_median = np.median(fwhm_s)
    for s in range(x_s.size):
        Rstar = fwhm_s[s]
        row,col = int(y_s[s]), int(x_s[s])
        if (flatwt_rc[row-5:row+6,col-5:col+6].sum() < 121) | (fwhm_s[s] <=0):
            Rstar = 1.5*fwhm_median                                    # if near an edge or bogus, use biggish star fwhm                                                               
        if fluxmax_s[s]/maxflatval > 0.4: Rstar *= -0.5*(np.log10(psflim/(fluxmax_s[s]/maxflatval)) + Apsf)/Bpsf
        flat_rc[np.where(np.sqrt((x_s[s]-1-x_im)**2+(y_s[s]-1-y_im)**2) < Rstar)] = 0
    os.remove(outtxt)
    os.remove(wtfile)        
#    pyfits.PrimaryHDU(flat_rc.astype("float32")).writeto(fits.split("/")[-1].split(".")[0] + "_sky.fits",clobber=True)

# separate amplifiers and correct for gain
    gain_a,ampcol_ac,maskcol_c = ampgains(flat_rc,rcbin_d)
    for amp in range(6): flat_rc[:,ampcol_ac[amp]] /= gain_a[amp]

# find edge of fov by finding peaks in 2nd derivative of row and column profiles, and mask beyond it
    edge_db = np.zeros((2,2))
    for d in (0,1):
        flat_b = flat_rc.sum(axis=1-d)                                                  # row, column sums
        flat_b[flat_b>0] /= (flat_rc>0).sum(axis=1-d)[flat_b>0]                         # r,c average of nonzero elements
        edge_db[d,0] = np.abs(flat_b[1:shape_d[d]/2]-flat_b.max()/10.).argmin() + 1
        edge_db[d,1] = np.abs(flat_b[shape_d[d]/2:-1]-flat_b.max()/10.).argmin() + shape_d[d]/2
    pixcen_d = edge_db.mean(axis=1)*rcbin_d[d]                                          # 1st guess: rc sum drops to 10% max

    pixrad = 240.0/pix_scale                                                            # 1st guess radius in pixels
    drad = -2./pix_scale                                                         # safety reduction of fov radius, in pixels  

    slice_d = np.empty(2,dtype=object)
    for d in (0,1):
        slice_d[d] = slice(shape_d[d])
        slice_d[1-d] = slice((pixcen_d[1-d]-32)/rcbin_d[1-d], (pixcen_d[1-d]+32)/rcbin_d[1-d])
        profile = tuple(slice_d)
        flat_b = flat_rc[profile].sum(axis=1-d)                                         # row, column profiles over center 64 pix
        flat_b[flat_b>0] /= (flat_rc>0)[profile].sum(axis=1-d)[flat_b>0]                # r,c average of nonzero elements
        okbin = np.zeros_like(flat_b,dtype=bool)
        okbin[1:-1] = (flat_b[:-2]>0) & (flat_b[1:-1]>0) & (flat_b[2:]>0)
        flat_b[1:-1] -= 0.5*(flat_b[0:-2] + flat_b[2:])                                 # 2nd derivative
        flat_b[~okbin] = 0.
        edge_db[d,0] = flat_b[max(0,(pixcen_d[d]-1.1*pixrad)/rcbin_d[d]):(pixcen_d[d]-0.9*pixrad)/rcbin_d[d]].argmax() \
                            + max(0,(pixcen_d[d]-1.1*pixrad)/rcbin_d[d])
        edge_db[d,1] = flat_b[(pixcen_d[d]+0.9*pixrad)/rcbin_d[d]:min(shape_d[d],(pixcen_d[d]+1.1*pixrad)/rcbin_d[d])].argmax() \
                            + (pixcen_d[d]+0.9*pixrad)/rcbin_d[d]
    pixcen_d = edge_db.mean(axis=1)*rcbin_d[d] 
    pixrad =  0.5*(edge_db[:,1]*rcbin_d[d] - edge_db[:,0]*rcbin_d[d]).mean() + drad
    flat_rc[((np.indices(shape_d)*rcbin_d[:,None,None] - pixcen_d[:,None,None])**2).sum(axis=0) > pixrad**2] = 0.   # zero out beyond fov edge
      
# find low-frequency shape by blocking to 64 unbinned pixels, then fitting with spline.  
    blk = 64/rcbin_d.min()
    flat_rc = blksmooth2d(flat_rc,blk,maskcol_c)
    flat_rc /= flat_rc[tuple(shape_d/2)]                                                # normalize to 1 at center                                          

# put gain back in
    if gaincor: 
        for amp in range(6): flat_rc[:,ampcol_ac[amp]] *= gain_a[amp]
    return flat_rc, pixcen_d

# ---------------------------------------------------------------------------------
def modelflat_rssimage(fits,pixcen_d): 
# Compute RSS imaging flatfield for a specific image
# _m = _yx tracker coordinates (m) (square array)
# _ae telescope focal plane azimuth,elevation coordinates (arcsec) (square array)
# _b = _rc instrument detector image row,column coordinates(bins),centered at pixcen_d

# adjustable parameters

    flddia = 540.       # diameter of flats in grid (arcsec)
    fldstop = 510.      # arcsec
    pixel = 15.             # pixel size in microns

    imgdist=np.loadtxt(distfile,usecols=(1,2))
    rc0_d = imgdist[0,::-1]
    rcd_d = imgdist[1,::-1]
    rd, A,B = imgdist[2,0],imgdist[3,0],imgdist[3,1]
    Fsclpoly=imgdist[4: ,0]

# input image and get geometry
    hdulist_image = pyfits.open(fits)
    image = fits.split(".")[-2][-4:]
    undist = fits.split(".")[-2].split("_")[-1] == "undist"       # *_undist.fits turns off distortion correction of model
    shape_d = np.asarray(hdulist_image["SCI"].data.shape)
    rcbin_d = np.array(hdulist_image[0].header["CCDSUM"].split(" ")).astype(int)
    rbin,cbin = rcbin_d
    dtype_image = hdulist_image['SCI'].data.dtype

    trkx = float(hdulist_image[0].header["TRKX"])
    trky = float(hdulist_image[0].header["TRKY"])
    trkrho = float(hdulist_image[0].header["TRKRHO"])

# input flat grid data for upper right track quadrant
# populate all quadrants by symmetry
    hdulist_model = pyfits.open(trackfile)
    flatgrid_yxea = hdulist_model[0].data
    eagrid = float(hdulist_model[0].header["CDELT1"])
    easize = flatgrid_yxea.shape[-1]
    yxgrid = float(hdulist_model[0].header["CDELT3"])
    yxsize = flatgrid_yxea.shape[0]
    flatgrid_yxea = np.concatenate((flatgrid_yxea[::-1,:,::-1,:],flatgrid_yxea[1:yxsize,:,:,:]),axis=0)
    flatgrid_yxea = np.concatenate((flatgrid_yxea[:,::-1,:,::-1],flatgrid_yxea[:,1:yxsize,:,:]),axis=1)
    trkpts_ea = (flatgrid_yxea>0).astype(int).sum(axis=0).sum(axis=0)
    maxtrkpts = trkpts_ea.max()
    flatgrid_yxea[:,:,(trkpts_ea < maxtrkpts)] = 0.
    oktrk = flatgrid_yxea.sum(axis=-1).sum(axis=-1) > 0.
    flatidx_hd = np.array(np.where(flatgrid_yxea.sum(axis=0).sum(axis=0) > 0.)).transpose()
    t_dt = yxgrid*np.arange(-yxsize+1,yxsize)
    flat_ea = np.zeros((easize,easize))

# interpolate in track space to image track location
    for (j,i) in flatidx_hd:
        trackint_h = ip.interp2d(t_dt,t_dt,flatgrid_yxea[:,:,j,i])
        flat_ea[j,i] = trackint_h(trkx,trky)

# prepare grid of location of RSS bin centers, and interpolate in field space, giving flat

    c = np.cos(np.radians(trkrho+180))                      # 180 is image inversion SALT -> RSS
    s = np.sin(np.radians(trkrho+180))
    rotate = np.transpose([[c, s],[-1.*s, c]])

    ji_drc = np.indices(shape_d)                            # output bin indices
    if undist:      
        Fsclpoly[5] = 120.
        ea_drc = np.indices(shape_d) - shape_d[:,None,None]/2   # bins, centered on optic axis
        ea_drc *= (pixel/Fsclpoly[5])*rcbin_d[:,None,None]      # RSS arcsec, centered on optic axis
        fcor_rc = np.ones_like(ea_drc[0])
    else:         
        rc_drc = np.indices(shape_d+1) - 0.5                    # start of bin corner grid
        rc_drc -= (rcd_d/rcbin_d+shape_d/2)[:,None,None]        # distorted bin corners, relative to distortion axis 
        pdist_rc = np.sqrt((rbin*rc_drc[0])**2 + (cbin*rc_drc[1])**2)    # distorted distance in pixels from distortion center
        p_rc = np.copy(pdist_rc)
        for iter in range(3):                                   # 3 iterations to undo distortion 
            dist_rc = (1. + A*(p_rc/rd)**2 + B*(p_rc/rd)**4)
            p_rc = pdist_rc/dist_rc
        rc_drc /= dist_rc[None,:,:]
        rc_drc += (rcd_d/rcbin_d)[:,None,None]                  # back to optic axis        
        ea_drc = (rc_drc[:,1:,1:] + rc_drc[:,1:,:-1] + rc_drc[:,:-1,1:] + rc_drc[:,:-1,:-1])/4.      # undistorted bin centers
        fcor_rc = (rc_drc[0,1:,:-1] - rc_drc[0,:-1,:-1]) * (rc_drc[1,:-1,1:] - rc_drc[1,:-1,:-1])    # approx bin size
        ea_drc *= (pixel/Fsclpoly[5])*rcbin_d[:,None,None]      # RSS arcsec (assume 5000 Ang, 7.5C)

    ea_drc = np.transpose(np.dot(np.transpose(ea_drc,(1,2,0)),rotate),(2,0,1))  # SALT arcsec
    infov_rc = np.sqrt(ea_drc[0]**2 + ea_drc[1]**2) < fldstop/2.
    ea_drc = (ea_drc + flddia/2.)/eagrid                        # bin corner fld grid coordinates
    
# update header later
    fov_dea = np.where(flat_ea>0)
    if fov_dea[0].shape[0] == 0:
        trkr = np.sqrt(trkx**2 + trky**2)
        print "Image %s Track Position X = %5.2f, Y = %5.2f, R = %5.2f  is beyond track flat model" % (image,trkx,trky,trkr)
        exit()
    flat_b = ip.griddata(fov_dea,flat_ea[fov_dea],tuple(ea_drc[:]),method="cubic",fill_value=0.)
    flat_rc = np.nan_to_num(flat_b.reshape(shape_d)*infov_rc*fcor_rc)
    flat_rc /= flat_rc[tuple(shape_d/2)]                          # normalize to 1 at center
    flat_rc =  img_recenter(flat_rc,pixcen_d/rcbin_d)             # recenter as requested

    return flat_rc

# ---------------------------------------------------------------------------------
def flatfld_rssimage(lampfits,skyfits,targlist,option=""):
    
    hdulist_image = pyfits.open(targlist[0])
    target_rc = np.copy(hdulist_image['SCI'].data)
    dtype_image = hdulist_image['SCI'].data.dtype
    rcbin_d = np.array(hdulist_image[0].header["CCDSUM"].split(" ")).astype(int)
    shape_d = np.array(target_rc.shape)
    pixcen_d = rcbin_d*shape_d/2.
    flat0_rc = np.ones_like(target_rc)

    print "RSS Imaging Flatfield using telescope model ",trkgrid

    if skyfits != "none":
        if skyfits.split(".")[-2].split("_")[-1] == "undist":
            print "Sky flats not supported for scrunched data"
            exit()  
        flat0_rc, pixcen_d = skyflat_rssimage(skyfits)
        skymodel_rc = modelflat_rssimage(skyfits,pixcen_d)
        modelpix = skymodel_rc > 0
    if lampfits != "none":
        if lampfits.split(".")[-2].split("_")[-1] == "undist":
            print "Lamp flats not supported for scrunched data"
            exit()  
        flat0_rc *= lampflat_rssimage(lampfits,gaincor=(skyfits=="none"))

    for targfits in targlist:
        hdulist_image = pyfits.open(targfits)
        target_rc = np.copy(hdulist_image["SCI"].data)
        flat_rc = flat0_rc*modelflat_rssimage(targfits,pixcen_d)
        if skyfits != "none":
            flat_rc[modelpix] /= skymodel_rc[modelpix]
#        flat_rc = flat0_rc                                          # test case with no track correction
        hdulist_image["SCI"].data = np.zeros_like(target_rc)
        flatpix = flat_rc > 0.
        hdulist_image["SCI"].data[flatpix] = target_rc[flatpix]/flat_rc[flatpix].astype(dtype_image)
        history = 'FLATFLD: lampflat= '+lampfits+' skyflat= '+skyfits+' trkgrid= '+trkgrid
        hdulist_image[0].header.add_history(history)
        flattenedfile = "f"+targfits.split("/")[-1]
        hdulist_image.writeto(flattenedfile, clobber = "True")    
        print "saved :",flattenedfile,                       
        if len(option)>0:
            hdulist_image["SCI"].data = flat_rc.astype(dtype_image)
            flatfile = targfits.split("/")[-1].split(".")[0] + "_flat.fits"
            print flatfile, 
            hdulist_image.writeto(flatfile, clobber = "True")
        print
    return

# ---------------------------------------------------------------------------------
if __name__=='__main__':

    lampfits=sys.argv[1]
    skyfits=sys.argv[2]
    targlist=sys.argv[3:]
    option = ""
    if (targlist[-1] == "saveflat") | (targlist[-1] == "debug"):
        option = targlist.pop()

    flatfld_rssimage(lampfits,skyfits,targlist,option)
