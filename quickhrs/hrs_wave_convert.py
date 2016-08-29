

import sys

from astropy.io import fits
from astropy import units as u
from astropy import constants as c
from PySpectrograph.Spectra.Spectrum import air2vac

from pyraf import iraf
from iraf import astutil

    

def convert_data(hdu, vslr):
   """Convert wavelenght array for vacumm wavlength and to the v_slr 
     
   """
   wave = hdu[1].data['Wavelength']
   wave = air2vac(wave)
   return wave * (1+vslr/c.c)

def rv_correct(hdu):
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
   img = sys.argv[1]
   hdu = fits.open(img)
   vlsr = rv_correct(hdu) * u.km/u.s
   wave = convert_data(hdu, vlsr)
   hdu[1].data['Wavelength'] = wave
   hdu.writeto('t'+img, clobber=True)
