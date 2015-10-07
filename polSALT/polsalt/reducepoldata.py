import os, sys, glob
reddir = '/d/carol/Synched/software/SALT/sandboxcopy/polSALT/polsalt/'
propid = 'NORDSIECK'

sys.path.append(reddir)
datadir = reddir+'data/'
import numpy as np
import pyfits

from imred import imred

from specpolwavmap import specpolwavmap
from specpolextract import specpolextract
from specpolrawstokes import specpolrawstokes
from specpolfinalstokes import specpolfinalstokes

obsdate = sys.argv[1]

os.chdir(obsdate)
if not os.path.isdir('sci'): os.mkdir('sci')
os.chdir('sci')

#basic image reductions
infile_list = glob.glob('../raw/P*fits')
imred(infile_list, './', datadir+'bpm_rss_11.fits', cleanup=True)

#basic polarimetric reductions
logfile='specpol'+obsdate+'.log'

#target and wavelength map
infile_list = sorted(glob.glob('m*fits'))
specpolwavmap(infile_list, propid, inter=True, logfile=logfile)

#background subtraction and extraction
infile_list = sorted(glob.glob('w*fits'))
specpolextract(infile_list, propid, logfile=logfile)

#raw stokes
infile_list = sorted(glob.glob('e*fits'))
specpolrawstokes(infile_list, propid, logfile=logfile)

#final stokes
infile_list = sorted(glob.glob('*_h*.fits'))
specpolfinalstokes(infile_list, propid, logfile=logfile)
