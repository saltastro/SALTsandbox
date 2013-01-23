
"""Data reduction script for AGN data.  This includes steps that are not yet included in the pipeline 
for best reduction of data

"""

import os, sys, glob

sys.path.append('/Users/crawford/agn/data/src/agnred')

from agnred import agnred

rawdir='../raw/'
prodir=os.path.curdir+'/'
imreduce=True 
specreduce=False
calfile='../../0506/p2/LTT4364.20120506_3.sens'
agnred(rawdir, prodir, imreduce, specreduce, calfile)
