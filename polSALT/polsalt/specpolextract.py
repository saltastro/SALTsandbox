
"""
specpolextract

Background subtract and extract spectropolarimetric spectrum data.  The end of the image-specific pipeline.
Ignore spectral curvature to start

"""

import os, sys, glob, shutil, inspect

import numpy as np
import pyfits
from scipy.interpolate import interp1d
import reddir
datadir = os.path.dirname(inspect.getfile(reddir))+"/data/"

from pyraf import iraf
from iraf import pysalt
from saltobslog import obslog
from saltsafelog import logging
debug = True

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

        hdu0 =  pyfits.open(outfile_list[0])       
        rows,cols = hdu0['SCI'].data.shape[1:3]
        rbin,cbin = np.array(obs_dict["CCDSUM"][0].split(" ")).astype(int)
        lam_c = hdu0['WAV'].data[0,rows/2]
        files = len(outfile_list)

    # sum spectra to find target and its background
        count = 0
        for i in range(files):
            if count == 0: image_orc = pyfits.open(outfile_list[i])['SCI'].data
            else: image_orc += pyfits.open(infile_list[i])['SCI'].data
            count += 1
        if count ==0:
            print 'No valid images'
            exit()
        image_orc /= count
        
    #   sum center 1/8 cols to get cross; 
    #   look for O,E max in cross within 200 unbinned pixels (+/-16 arcsec) of expected posn
    #   take background +/- 112 - 128 unbinned pixels (+/-14-16 arcsec) from max
    #   find edge of spectrum where cross signal drops to 2%, edge of background where it deviates by more than 2%
    #   do target identification in transposed cr space to have boolean indexing flatten correctly

        lam_m = np.loadtxt(datadir+"wollaston.txt",dtype=float,usecols=(0,))
        rpix_om = np.loadtxt(datadir+"wollaston.txt",dtype=float,unpack=True,usecols=(1,2))
        cross_or = np.sum(image_orc[:,:,cols/2-cols/16:cols/2+cols/16],axis=2)
        image_ocr = np.transpose(image_orc,(0,2,1)) 

        srow_oc = np.zeros((2,cols))
        tgtspec_ocr = np.zeros((2,cols,rows),dtype=bool)
        tgtbkg_ocr = np.zeros((2,cols,rows),dtype=bool)
        col_cr,row_cr = np.indices(image_ocr.shape[1:3])

        for o in (0,1):
            expectrow = (1-o)*rows + interp1d(lam_m,rpix_om[o],kind='cubic')(lam_c[cols/2])/rbin    
            crossmaxval = np.max(cross_or[o,expectrow-100/rbin:expectrow+100/rbin])
            drow = np.where(cross_or[o]==crossmaxval)[0][0] - expectrow
            row_c = ((1-o)*rows + interp1d(lam_m,rpix_om[o],kind='cubic')(lam_c)/rbin + drow).astype(int)
            maxval_c =  np.max(image_ocr[o,np.abs(row_cr-row_c[:,None])<(20/rbin)].reshape((cols,-1)),axis=1)
            bkg1_c = np.median(image_ocr[o,np.abs(row_cr-row_c[:,None]-120)< 8/rbin].reshape((cols,-1)),axis=1)
            bkg2_c = np.median(image_ocr[o,np.abs(row_cr-row_c[:,None]+120)< 8/rbin].reshape((cols,-1)),axis=1)
            bkg_c = 0.5*(bkg1_c+bkg2_c)
            edge_c = 0.02*(maxval_c-bkg_c)+bkg_c
            imagebot_cr = image_ocr[o,(row_cr<row_c[:,None])&(row_cr>(row_c[:,None]-row_c.min()))].reshape((cols,-1))
            imagetop_cr = image_ocr[o,(row_cr>row_c[:,None])&(row_cr<(row_c[:,None]+(rows-row_c.max())))].reshape((cols,-1))
            esrow0_c = row_c - np.argmax(np.fliplr(imagebot_cr) < edge_c[:,None],axis=1)  
            esrow1_c = row_c + np.argmax(imagetop_cr < edge_c[:,None],axis=1)

            imagebot_cr = image_ocr[o,(row_cr<esrow0_c[:,None]-2)    \
                        &(row_cr>(esrow0_c[:,None]-esrow0_c.min()))].reshape((cols,-1))
            imagetop_cr = image_ocr[o,(row_cr>esrow1_c[:,None]+2)    \
                        &(row_cr<(esrow1_c[:,None]+(rows-esrow1_c.max())))].reshape((cols,-1))
            ebrow0_c = esrow0_c-3 - np.argmax(np.fliplr(imagebot_cr) > edge_c[:,None],axis=1)
            ebrow1_c = esrow1_c+3 + np.argmax(imagetop_cr > edge_c[:,None],axis=1)

            tgtspec_ocr[o] = (row_cr >= esrow0_c[:,None]) & (row_cr <= esrow1_c[:,None])
            tgtbkg_ocr[o] = ((row_cr >= ebrow0_c[:,None]) & (row_cr <  esrow0_c[:,None])) \
                |   ((row_cr >  esrow1_c[:,None]) & (row_cr <= ebrow1_c[:,None]))

            srow_oc[o] = (esrow0_c+esrow1_c)/2

        log.message('Target center row: O '+("%4i " % srow_oc[0,cols/2])    \
                    +'  E '+("%4i " % srow_oc[1,cols/2]), with_header=False)

        tgtspec_orc = np.transpose(tgtspec_ocr,(0,2,1))
        tgtbkg_orc = np.transpose(tgtbkg_ocr,(0,2,1))

    # extract background-subtracted spectra

        for i in range(files):
            hdulist = pyfits.open(outfile_list[i])
            sci_orc = hdulist['sci'].data
            var_orc = hdulist['var'].data
            bpm_orc = hdulist['bpm'].data
            wav_orc = hdulist['wav'].data 

        # output column flagged bad if bkg or target has no good data or target contains bad pix  
            specpix = tgtspec_orc & (bpm_orc==0)        
            bkgpix = tgtbkg_orc & (bpm_orc==0)
            bpm_oc = (bkgpix.sum(axis=1)==0) | (specpix.sum(axis=1)==0) \
                    | ((tgtspec_orc & (bpm_orc==1)).sum(axis=1) >0)
            bkg_oc = np.zeros((2,cols)); varbkg_oc = np.zeros_like(bkg_oc)
            bkg_oc[~bpm_oc] = (sci_orc*bkgpix).sum(axis=1)[~bpm_oc]/bkgpix.sum(axis=1)[~bpm_oc]
            varbkg_oc[~bpm_oc] = (var_orc*bkgpix).sum(axis=1)[~bpm_oc]/(bkgpix.sum(axis=1)[~bpm_oc])**2
            sci_oc = ((sci_orc-bkg_oc[:,None,:])*specpix).sum(axis=1)
            var_oc = ((var_orc-bkg_oc[:,None,:])*specpix).sum(axis=1) + specpix.sum(axis=1)*varbkg_oc
            tgt_oc = tgtspec_orc.sum(axis=1) > 0
            wav_oc = np.zeros_like(sci_oc)
            wav_oc[tgt_oc] = (wav_orc*tgtspec_orc).sum(axis=1)[tgt_oc] / tgtspec_orc.sum(axis=1)[tgt_oc]

        # write O,E spectrum, prefix "s". VAR, BPM, WAV for each spectrum. y dim is virtual (length 1)
        # for consistency with other modes

            hdulist['SCI'].data = sci_oc.reshape((2,1,cols)).astype('float32')
            hdulist['VAR'].data = var_oc.reshape((2,1,cols)).astype('float32')
            hdulist['BPM'].data = bpm_oc.reshape((2,1,cols)).astype('uint8')
            hdulist['WAV'].data = wav_oc.reshape((2,1,cols)).astype('float32')
            hdulist.writeto('e'+outfile_list[i],clobber=True,output_verify='warn')
            log.message('Output file '+'e'+outfile_list[i] , with_header=False)

    return

