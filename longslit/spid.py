#/usr/bin/env python
import sys
from pyraf import iraf
from pyraf.iraf import pysalt
from specidentify import specidentify

specidentify(images=sys.argv[1], linelist=sys.argv[2], outfile='db.sol', guesstype='rss', guessfile='', automethod='Matchlines', function='polynomial', order=3,rstep=100 , rstart='middlerow', mdiff=20, thresh=3,niter=5,inter=True, clobber=False, logfile='salt.log', verbose=True)
