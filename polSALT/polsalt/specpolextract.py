
"""
specpolextract

Background subtract and extract spectropolarimetric data.  The end of the image-specific pipeline.

"""

import os, sys, glob, shutil, inspect

import numpy as np
import pyfits
from scipy.interpolate import interp1d
from scipy import linalg as la

import reddir
datadir = os.path.dirname(inspect.getfile(reddir))+"/data/"

from blksmooth2d import blksmooth2d
from specpollampextract import specpollampextract
from specpolsignalmap import specpolsignalmap
from skysub2d_khn import make_2d_skyspectrum
from scrunch1d import scrunch1d
from pyraf import iraf
from iraf import pysalt
from saltobslog import obslog
from saltsafelog import logging

# np.seterr(invalid='raise')
np.set_printoptions(threshold=np.nan)
debug = True

# ---------------------------------------------------------------------------------
def specpolextract(infile_list, propcode=None, logfile='salt.log'):

#set up the files
    obsdate=os.path.basename(infile_list[0])[8:16]

    with logging(logfile, debug) as log:
        #create the observation log
        obs_dict=obslog(infile_list)
        outfile_list = []

        for i in range(len(infile_list)):
            if (obs_dict['OBJECT'][i].upper().strip()=='ARC'): continue
            outfile_list.append(infile_list[i])

        obs_dict=obslog(outfile_list)
        hdu0 =  pyfits.open(outfile_list[0])       
        rows,cols = hdu0['SCI'].data.shape[1:3]
        rbin,cbin = np.array(obs_dict["CCDSUM"][0].split(" ")).astype(int)

# special version for lamp data
        lampid = obs_dict["LAMPID"][0].strip().upper()
        if lampid!="NONE":
            specpollampextract(outfile_list, logfile=logfile)           
            return

        slitid = obs_dict["MASKID"][0]
        if slitid[0] =="P": slitwidth = float(slitid[2:5])/10.
        else: slitwidth = float(slitid)

# sum spectra to find target, background artifacts, and estimate sky flat and psf functions
        count = 0
        files = len(outfile_list)
        for i in range(files):
            isbadbin_orc = pyfits.open(outfile_list[i])['BPM'].data > 0
            if count == 0: 
                count_orc = (~isbadbin_orc).astype(int)
                image_orc = pyfits.open(outfile_list[i])['SCI'].data*count_orc
                var_orc = pyfits.open(outfile_list[i])['VAR'].data
            else:
                count_orc += (~isbadbin_orc).astype(int)
                image_orc += pyfits.open(infile_list[i])['SCI'].data*(~isbadbin_orc).astype(int)
                var_orc += pyfits.open(outfile_list[i])['VAR'].data
            count += 1
        if count ==0:
            print 'No valid images'
            exit()
        image_orc[count_orc>0] /= count_orc[count_orc>0]
        isbadbin_orc = (count_orc==0) | (image_orc==0)
        isbadbinpol_orc = (count_orc < count) | (image_orc==0)    # conservative bpm for pol extraction
        var_orc[count_orc>0] /= count_orc[count_orc>0]**2
        wav_orc = pyfits.open(outfile_list[0])['WAV'].data 

        hdusum = pyfits.PrimaryHDU(header=hdu0[0].header)    
        hdusum = pyfits.HDUList(hdusum)
        header=hdu0['SCI'].header.copy()           
        hdusum.append(pyfits.ImageHDU(data=image_orc, header=header, name='SCI'))
        hdusum.append(pyfits.ImageHDU(data=var_orc, header=header, name='VAR'))
        hdusum.append(pyfits.ImageHDU(data=isbadbin_orc.astype('uint8'), header=header, name='BPM'))
        hdusum.append(pyfits.ImageHDU(data=wav_orc, header=header, name='WAV'))

        psf_orc,skyflat_orc,isnewbadbin_orc,isbkgcont_orc,maprow_od,drow_oc = \
                specpolsignalmap(hdusum,logfile=logfile)
  
        maprow_ocd = maprow_od[:,None,:] - drow_oc[:,:,None]

        isedge_orc = (np.arange(rows)[:,None] < maprow_ocd[:,None,:,0]) | \
            (np.arange(rows)[:,None] > maprow_ocd[:,None,:,3])
        istarget_orc = (np.arange(rows)[:,None] > maprow_ocd[:,None,:,1]) & \
            (np.arange(rows)[:,None] < maprow_ocd[:,None,:,2])
        isskycont_orc = (((np.arange(rows)[:,None] < maprow_ocd[:,None,:,0]+rows/16) |  \
            (np.arange(rows)[:,None] > maprow_ocd[:,None,:,3]-rows/16)) & ~isedge_orc)                                     
        isbkgcont_orc &= (~isbadbin_orc & ~isedge_orc & ~istarget_orc)
        isbadbin_orc |= isnewbadbin_orc
        isbadbinpol_orc |= isnewbadbin_orc

#        pyfits.PrimaryHDU(isbadbinpol_orc.astype('uint8')).writeto('isbadbinpol_orc.fits',clobber=True)  
#        pyfits.PrimaryHDU(isnewbadbin_orc.astype('uint8')).writeto('isnewbadbin_orc.fits',clobber=True)  

    # scrunch and normalize psf for optimized extraction
        psfnormmin = 0.70    # wavelengths with less than this flux in good bins are marked bad
        wbin = wav_orc[0,rows/2,cols/2]-wav_orc[0,rows/2,cols/2-1] 
        wbin = float(int(wbin/0.75))
        wmin,wmax = wav_orc.min(axis=2).max(),wav_orc.max(axis=2).min()
        wedgemin = wbin*int(wmin/wbin+0.5) + wbin/2.
        wedgemax = wbin*int(wmax/wbin-0.5) + wbin/2.
        wedge_w = np.arange(wedgemin,wedgemax+wbin,wbin)
        wavs = wedge_w.shape[0] - 1
        binedge_orw = np.zeros((2,rows,wavs+1))
        psf_orw = np.zeros((2,rows,wavs));  isbadbin_orw = np.zeros((2,rows,wavs),dtype=bool); 
        specrow_or = maprow_od[:,1:3].mean(axis=1)[:,None] + np.arange(-rows/4,rows/4)
        for o in (0,1):
            for r in specrow_or[o]:
                binedge_orw[o,r] = interp1d(wav_orc[o,r],np.arange(cols))(wedge_w)
                psf_orw[o,r] = scrunch1d(psf_orc[o,r],binedge_orw[o,r])
                isbadbin_orw[o,r] = scrunch1d(isbadbinpol_orc[o,r].astype(int),binedge_orw[o,r]) > 0.01
        okbin_orw = (isbadbin_orw == 0)
        psf_orw /= psf_orw.sum(axis=1)[:,None,:]
        okbin_orw &= ((psf_orw*okbin_orw).sum(axis=1)[:,None,:] > psfnormmin)

#        np.savetxt("psfnorm_ow.txt",(psf_orw*okbin_orw).sum(axis=1).T,fmt="%10.4f") 
#        pyfits.PrimaryHDU(okbin_orw.astype('uint8')).writeto('okbin_orw.fits',clobber=True)  
#        pyfits.PrimaryHDU(psf_orw.astype('float32')).writeto('psf_orw.fits',clobber=True)    

    # background-subtract and extract spectra

        for i in range(files):
            hdulist = pyfits.open(outfile_list[i])
            sci_orc = hdulist['sci'].data
            var_orc = hdulist['var'].data
            tnum = os.path.basename(outfile_list[i]).split('.')[0][-3:]

        # make background continuum image, smoothed over resolution element
            rblk,cblk = int(1.5*8./rbin), int(slitwidth*8./cbin)
            target_orc = np.zeros_like(sci_orc)    
            for o in (0,1):
                bkgcont_rc = blksmooth2d(sci_orc[o],isbkgcont_orc[o],rblk,cblk,0.25,mode="mean")           
                    
        # remove sky continuum: ends of bkg continuum * skyflat
                skycont_c = (bkgcont_rc.T[isskycont_orc[o].T]/skyflat_orc[o].T[isskycont_orc[o].T])  \
                            .reshape((cols,-1)).mean(axis=1)
                skycont_rc = skycont_c*skyflat_orc[o]
        
        # remove sky lines: image - bkg cont run through 2d sky averaging
                obj_data = ((sci_orc[o] - bkgcont_rc)/skyflat_orc)[o]
                obj_data[(isbadbin_orc | isedge_orc | istarget_orc)[o]] = np.nan
#                pyfits.PrimaryHDU(obj_data.astype('float32')).writeto('obj_data.fits',clobber=True)
                skylines_rc = make_2d_skyspectrum(obj_data,wav_orc[o],np.array([[0,rows],]))*skyflat_orc[o]
                target_orc[o] = sci_orc[o] - skycont_rc - skylines_rc
#                pyfits.PrimaryHDU(skylines_rc.astype('float32')).writeto('skylines_rc_'+tnum+'_'+str(o)+'.fits',clobber=True)
            target_orc *= (~isbadbinpol_orc).astype(int)             
#            pyfits.PrimaryHDU(target_orc.astype('float32')).writeto('target_'+tnum+'_orc.fits',clobber=True)

        # extract spectrum optimally (Horne, PASP 1986)
            target_orw = np.zeros((2,rows,wavs));   var_orw = np.zeros_like(target_orw)
            for o in (0,1):
                for r in specrow_or[o]:
                    target_orw[o,r] = scrunch1d(target_orc[o,r],binedge_orw[o,r])
                    var_orw[o,r] = scrunch1d(var_orc[o,r],binedge_orw[o,r])
#            pyfits.PrimaryHDU(target_orw.astype('float32')).writeto('target_'+tnum+'_orw.fits',clobber=True)
            ok_orw = (okbin_orw &(var_orw > 0))

        # use master psf shifted in row to allow for guide errors
            pwidth = 2*int(1./psf_orw.max())
            ok_w = ((psf_orw*(~ok_orw)).sum(axis=1) < 0.03/float(pwidth/2)).all(axis=0)
            crosscor_s = np.zeros(pwidth)
            for s in range(pwidth):
                crosscor_s[s] = (psf_orw[:,s:s-pwidth]*target_orw[:,pwidth/2:-pwidth/2]*ok_w).sum()
            smax = np.argmax(crosscor_s)
            s_S = np.arange(smax-pwidth/4,smax-pwidth/4+pwidth/2+1)
            polycof = la.lstsq(np.vstack((s_S**2,s_S,np.ones_like(s_S))).T,crosscor_s[s_S])[0]
            pshift = -(-0.5*polycof[1]/polycof[0] - pwidth/2)
            s = int(pshift+pwidth)-pwidth
            sfrac = pshift-s
            psfsh_orw = np.zeros_like(psf_orw)
            outrow = np.arange(max(0,s+1),rows-(1+int(abs(pshift)))+max(0,s+1))

            psfsh_orw[:,outrow] = (1.-sfrac)*psf_orw[:,outrow-s] + sfrac*psf_orw[:,outrow-s-1]  
            wt_orw = np.zeros_like(target_orw)
            wt_orw[ok_orw] = psfsh_orw[ok_orw]/var_orw[ok_orw]
            var_ow = (psfsh_orw*wt_orw).sum(axis=1)
            ok_ow = (var_ow > 0)
            bpm_ow = (~ok_ow).astype('uint8')
            var_ow[ok_ow] = 1./var_ow[ok_ow]
            sci_ow = var_ow*(target_orw*wt_orw).sum(axis=1)

        # write O,E spectrum, prefix "s". VAR, BPM for each spectrum. y dim is virtual (length 1)
        # for consistency with other modes
            hduout = pyfits.PrimaryHDU(header=hdulist[0].header)    
            hduout = pyfits.HDUList(hduout)
            header=hdulist['SCI'].header.copy()
            header.update('VAREXT',2)
            header.update('BPMEXT',3)
            header.update('CRVAL1',wedge_w[0]+wbin/2.)
            header.update('CRVAL2',0)
            header.update('CDELT1',wbin)
            header.update('CTYPE1','Angstroms')
            
            hduout.append(pyfits.ImageHDU(data=sci_ow.reshape((2,1,wavs)), header=header, name='SCI'))
            header.update('SCIEXT',1,'Extension for Science Frame',before='VAREXT')
            hduout.append(pyfits.ImageHDU(data=var_ow.reshape((2,1,wavs)), header=header, name='VAR'))
            hduout.append(pyfits.ImageHDU(data=bpm_ow.reshape((2,1,wavs)), header=header, name='BPM'))            
            
            hduout.writeto('e'+outfile_list[i],clobber=True,output_verify='warn')
            log.message('Output file '+'e'+outfile_list[i] , with_header=False)
      
    return

