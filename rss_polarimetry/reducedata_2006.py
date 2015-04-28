import os, sys, glob
reddir = '[Update with path to zSALT]/zSALT/zsalt/'
sys.path.append(reddir)
bpmfile = '[update wtih path to bpm]/bpm.fits'
propid = '0'

from imred import imred
from specred import specred

from astropy.io import fits


ddir = sys.argv[1]

os.chdir(ddir)
if not os.path.isdir('sci'): os.mkdir('sci')
os.chdir('sci')

#basic image reuctions
infile_list = glob.glob('../raw/P*fits')

for img in infile_list:
    hdu = fits.open(img, 'update')
    hdu[1].header['XTALK'] = 1474
    hdu[2].header['XTALK'] = 1474
    hdu[3].header['XTALK'] = 1166
    hdu[4].header['XTALK'] = 1111
    hdu[5].header['XTALK'] = 1377
    hdu[6].header['XTALK'] = 1377
    hdu.close()


imred(infile_list, './', bpmfile, cleanup=True)

#spectroscopic reductions
infile_list = glob.glob('m*fits')
specred(infile_list, propid, inter=True)

#UNEDIT IF YOU WANT TO REDUCE MOS DATA
#mosxml = [MOSXMLFILE]
#mosred = infile_list, propid, mosxml, inter=True)





