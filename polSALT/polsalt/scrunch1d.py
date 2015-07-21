#! /usr/bin/env python

# Resample data into new bins, preserving flux

import os, sys, time, glob, shutil
import numpy as np
import pyfits
WCSHdr = ['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2','CTYPE1','CTYPE2','CUNIT1','CUNIT2']

def scrunch1d(input_a,binedge_x):

    na = input_a.size
    nx = binedge_x.size - 1

    inputT_a = np.append(input_a,0)                         # deal with edge of array
    binedgeT_x = np.maximum(0.,np.minimum(na,binedge_x))    # deal with edge of array

    ns = na+nx+1
    ia_s = np.append(binedgeT_x,range(na+1))
    ix_s = np.zeros_like(ia_s) - 1
    ix_s[0:nx+1] = np.arange(nx+1)
    ia_argsort_s = np.argsort(ia_s)
    ia_s = ia_s[ia_argsort_s]
    origin_s = ia_s.astype(int)[0:ns]
    value_s = inputT_a[origin_s]*(ia_s[1:ns+1] - ia_s[0:ns])
    destination_s = ix_s[ia_argsort_s].astype(int)[0:ns]
    destination_s = np.maximum.accumulate(destination_s)
    ix_x = np.arange(nx)    
    output_x = np.where(destination_s == ix_x[:,None],value_s,0).sum(axis=1)
    return output_x

def scrunch1dtest(inputfile,binedgefile):
    input_a = np.loadtxt(inputfile)
    print input_a, input_a.sum()
    binedge_x = np.loadtxt(binedgefile)
    print binedge_x
    output_x = scrunch1D(input_a,binedge_x)
    print output_x, output_x.sum()
    return

if __name__=='__main__':
    inputfile=sys.argv[1]
    binedgefile=sys.argv[2]
    scrunch1dtest(inputfile,binedgefile)

