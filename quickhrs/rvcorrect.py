

import sys
from astropy.io import fits
    
from pyraf import iraf
from iraf import astutil



def rv_correct(img):
    hdu = fits.open(img)
    date_obs = hdu[0].header['DATE-OBS'].split('-')
    iraf.rvcorrect(header="N", input="N",
                   imupdate="N",
                   observatory="SAAO",
                   year=float(date_obs[0]),
                   month=float(date_obs[1]),
                   day=float(date_obs[2]),
                   ut=hdu[0].header['TIME-OBS'],
                   ra=hdu[0].header['RA'],
                   dec=hdu[0].header['DEC'],
                   vobs=0,
                   mode="a")
    return iraf.rvcorrect.vlsr

if __name__=='__main__':
   rv_correct(sys.argv[1])
   
