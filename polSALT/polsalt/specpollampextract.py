
"""
specpollampextract

Extract spectropolarimetric lamp spectrum data.  

"""

import os, sys, glob, shutil, inspect

import numpy as np
import pyfits
from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import shift

import reddir
datadir = os.path.dirname(inspect.getfile(reddir))+"/data/"

from scrunch1d import scrunch1d
from pyraf import iraf
from iraf import pysalt
from saltobslog import obslog
from saltsafelog import logging

# np.seterr(invalid='raise')
np.set_printoptions(threshold=np.nan)
debug = True

# ---------------------------------------------------------------------------------
def specpollampextract(infile_list, logfile='salt.log'):

    obsdate=os.path.basename(infile_list[0])[8:16]

    with logging(logfile, debug) as log:
        log.message('Extraction of Lamp Images' , with_header=False)
        obs_dict=obslog(infile_list)

        hdu0 =  pyfits.open(infile_list[0])       
        rows,cols = hdu0['SCI'].data.shape[1:3]
        rbin,cbin = np.array(obs_dict["CCDSUM"][0].split(" ")).astype(int)
        slitid = obs_dict["MASKID"][0]

        lam_c = hdu0['WAV'].data[0,rows/2]
        files = len(infile_list)
        outfile_list = infile_list

# sum spectra to find target
        count = 0
        for i in range(files):
            badbin_orc = pyfits.open(outfile_list[i])['BPM'].data.astype(bool)
            if count == 0: 
                count_orc = (~badbin_orc).astype(int)
                image_orc = pyfits.open(outfile_list[i])['SCI'].data*count_orc
                var_orc = pyfits.open(outfile_list[i])['VAR'].data
            else:
                count_orc += (~badbin_orc).astype(int)
                image_orc += pyfits.open(outfile_list[i])['SCI'].data*(~badbin_orc)
                var_orc += pyfits.open(outfile_list[i])['VAR'].data
            count += 1
        if count ==0:
            print 'No valid images'
            exit()
        image_orc[count_orc>0] /= count_orc[count_orc>0]
        badbin_orc = (count_orc==0) | (image_orc==0)
        badbinpol_orc = (count_orc < count) | (image_orc==0)    # conservative bpm for pol extraction
        var_orc[count_orc>0] /= count_orc[count_orc>0]**2
        wav_orc = pyfits.open(outfile_list[0])['WAV'].data 

        lam_m = np.loadtxt(datadir+"wollaston.txt",dtype=float,usecols=(0,))
        rpix_om = np.loadtxt(datadir+"wollaston.txt",dtype=float,unpack=True,usecols=(1,2))

    # trace spectrum, compute spatial profile 
        profile_orc = np.zeros_like(image_orc)
        drow_oc = np.zeros((2,cols))
        expectrow_oc = np.zeros((2,cols),dtype='float32')
        maxrow_oc = np.zeros((2,cols),dtype=int); maxval_oc = np.zeros((2,cols),dtype='float32')
        col_cr,row_cr = np.indices(image_orc[0].T.shape)
        cross_or = np.sum(image_orc[:,:,cols/2-cols/16:cols/2+cols/16],axis=2)
        
        okprof_oyc = np.ones((2,rows,cols),dtype='bool')
        profile_oyc = np.zeros_like(profile_orc)
        badbin_oyc = np.zeros_like(badbin_orc)
        var_oyc = np.zeros_like(var_orc)
        wav_oyc = np.zeros_like(wav_orc)

        for o in (0,1):
        # find spectrum roughly from max of central cut, then within narrow curved aperture
            expectrow_oc[o] = (1-o)*rows + interp1d(lam_m,rpix_om[o],kind='cubic')(lam_c)/rbin    
            crossmaxval = np.max(cross_or[o,expectrow_oc[o,cols/2]-100/rbin:expectrow_oc[o,cols/2]+100/rbin])
            drow = np.where(cross_or[o]==crossmaxval)[0][0] - expectrow_oc[o,cols/2]
            row_c = (expectrow_oc[o] + drow).astype(int)
            aperture_cr = ((row_cr-row_c[:,None])>=-20/rbin) & ((row_cr-row_c[:,None])<=20/rbin)            
            maxrow_oc[o] = np.argmax(image_orc[o].T[aperture_cr].reshape((cols,-1)),axis=1) + row_c - 20/rbin
            maxval_oc[o] = image_orc[o,maxrow_oc[o]].diagonal()
            drow1_c = maxrow_oc[o] -(expectrow_oc[o] + drow)
            trow_o = maxrow_oc[:,cols/2]

        # divide out spectrum (allowing for spectral curvature) to make spatial profile
            okprof_c = (maxval_oc[o] != 0)
            drow2_c = np.polyval(np.polyfit(np.where(okprof_c)[0],drow1_c[okprof_c],3),(range(cols)))
            okprof_c &= np.abs(drow2_c - drow1_c) < 3
            norm_rc = np.zeros((rows,cols))
            for r in range(rows):
                norm_rc[r] = interp1d(wav_orc[o,trow_o[o],okprof_c],maxval_oc[o,okprof_c], \
                    bounds_error = False, fill_value=0.)(wav_orc[o,r])
            okprof_rc = (norm_rc != 0.)
            profile_orc[o,okprof_rc] = image_orc[o,okprof_rc]/norm_rc[okprof_rc]
            var_orc[o,okprof_rc] = var_orc[o,okprof_rc]/norm_rc[okprof_rc]**2
            drow_oc[o] = -(expectrow_oc[o] - expectrow_oc[o,cols/2] + drow2_c -drow2_c[cols/2])      
            
    # Find, mask off fov edge, including possible beam overlap
        profile_or = np.median(profile_orc[:,:,cols/2-32:cols/2+32],axis=2)
        edgerow_do = np.zeros((2,2),dtype=int)
        badrow_or = np.zeros((2,rows),dtype=bool)
        axisrow_o = np.zeros(2)
#        np.savetxt("profile_or.txt",profile_or.T,fmt="%10.6f")

        maxoverlaprows = 34/rbin                       
        for d,o in np.ndindex(2,2):                         # _d = (0,1) = (bottom,top)
            row_y = np.where((d==1) ^ (np.arange(rows) < trow_o[o]))[0][::2*d-1]
            edgeval = np.median(profile_or[o,row_y],axis=-1)        
            hist,bin = np.histogram(profile_or[o,row_y],bins=32,range=(0,edgeval))
            histarg = 32 - np.argmax(hist[::-1]<3)      # edge: <3 in hist in decreasing dirn
            edgeval = bin[histarg]
            edgerow_do[d,o] = trow_o[o] + (2*d-1)*(np.argmax(profile_or[o,row_y] <= edgeval))
            axisrow_o[o] += edgerow_do[d,o]
            edgerow_do[d,o] = np.clip(edgerow_do[d,o],maxoverlaprows,rows-maxoverlaprows)
            badrow_or[o] |= ((d==1) ^ (np.arange(rows) < (edgerow_do[d,o]+d)))       
        axisrow_o /= 2.

    # find edge of target slit, if multiple slits
        slitrow_os = 0.5*(np.argwhere((profile_or[:,:-1] > 0.1) & (profile_or[:,1:] < 0.1)) + \
            np.argwhere((profile_or[:,:-1] < 0.1) & (profile_or[:,1:] > 0.1)))[:,1].reshape(2,-1)
        slits = slitrow_os.shape[1]
        targetslit = np.where(abs(trow_o[0] - slitrow_os[0]) < 4)[0][0]
        for o in (0,1):
            if targetslit > 0:
                edgerow_do[0,o] = slitrow_os[o,targetslit-1:targetslit+1].mean()
            if targetslit < slits-1:
                edgerow_do[1,o] = slitrow_os[o,targetslit:targetslit+2].mean()
        edgerow_doc = (edgerow_do[:,:,None] + drow_oc[None,:,:]).astype(int)
        bkgrows_do = ((trow_o - edgerow_do)/2).astype(int)
        bkgrow_doc = edgerow_doc + bkgrows_do[:,:,None]/2
        isbkg_dorc = (((np.arange(rows)[:,None] - edgerow_doc[:,:,None,:]) * \
              (np.arange(rows)[:,None] - edgerow_doc[:,:,None,:] + bkgrows_do[:,:,None,None])) < 0)
        if targetslit == 0:        
            log.message('Optical axis row:  O    %4i     E    %4i' % tuple(axisrow_o), with_header=False)
        log.message('Target center row: O    %4i     E    %4i' % tuple(trow_o), with_header=False)
        log.message('Bottom, top row:   O %4i %4i   E %4i %4i \n' \
                % tuple(edgerow_do.T.flatten()), with_header=False)
  
    # background-subtract and extract spectra

        for i in range(files):
            hdulist = pyfits.open(outfile_list[i])
            sci_orc = hdulist['sci'].data
            var_orc = hdulist['var'].data

        # make background continuum image, linearly interpolated in row
            bkg_doc = np.zeros((2,2,cols))
            for d,o in np.ndindex(2,2):
                bkg_doc[d,o] = np.median(sci_orc[o].T[isbkg_dorc[d,o].T].reshape((cols,-1)),axis=1)         
            bkgslp_oc = (bkg_doc[1] - bkg_doc[0])/(bkgrow_doc[1] - bkgrow_doc[0])
            bkgbase_oc = (bkg_doc[1] + bkg_doc[0])/2. - bkgslp_oc*(bkgrow_doc[1] + bkgrow_doc[0])/2.
            np.savetxt('bkg.txt',np.vstack((bkg_doc.reshape((4,-1)),bkgslp_oc,bkgbase_oc)).T,fmt="%11.4f")
            bkg_orc = bkgbase_oc[:,None,:] + bkgslp_oc[:,None,:]*np.arange(rows)[:,None]                  
            target_orc = sci_orc-bkg_orc             

        # extract spectrum 
            wbin = wav_orc[0,rows/2,cols/2]-wav_orc[0,rows/2,cols/2-1] 
            wbin = float(int(wbin/0.75))
            wmin,wmax = wav_orc.min(axis=2).max(),wav_orc.max(axis=2).min()
            wedgemin = wbin*int(wmin/wbin+0.5) + wbin/2.
            wedgemax = wbin*int(wmax/wbin-0.5) + wbin/2.
            wedge_w = np.arange(wedgemin,wedgemax+wbin,wbin)
            wavs = wedge_w.shape[0] - 1
            binedge_orw = np.zeros((2,rows,wavs+1))
            target_orw = np.zeros((2,rows,wavs));   var_orw = np.zeros_like(target_orw)
            badbin_orw = np.ones((2,rows,wavs))
            for o in (0,1):
                for r in range(edgerow_doc[0,o].min(),edgerow_doc[1,o].max()):
                    binedge_orw[o,r] = interp1d(wav_orc[o,r],np.arange(cols))(wedge_w)
                    target_orw[o,r] = scrunch1d(target_orc[o,r],binedge_orw[o,r])
                    var_orw[o,r] = scrunch1d(var_orc[o,r],binedge_orw[o,r])
                    badbin_orw[o,r] = scrunch1d(badbinpol_orc[o,r],binedge_orw[o,r])
            okbin_orw = (badbin_orw==0.) & \
                ((badbin_orw==0.).sum(axis=1) > (edgerow_do[1]-edgerow_do[0])[:,None]-3)[:,None,:]
#            pyfits.PrimaryHDU(badbin_orw.astype('float32')).writeto('badbin_orw.fits',clobber=True) 
#            pyfits.PrimaryHDU(okbin_orw.astype('uint8')).writeto('okbin_orw.fits',clobber=True) 

            sci_ow = (okbin_orw.astype(int)*target_orw).sum(axis=1).astype('float32')
            var_ow = (okbin_orw.astype(int)*var_orw).sum(axis=1).astype('float32')
            bpm_ow = (~okbin_orw).all(axis=1).astype('uint8')

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

