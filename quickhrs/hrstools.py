
from pylab import *

from scipy import ndimage as nd

from astropy.io import fits
from astropy import modeling as md

from hrsmodel import HRSModel

red_edge=5500

def get_onlyobject(data, obj=2, dy=20):
   """Return an array of the same shape as data
      but with only the object data and all 
      other data masked
   """
   ys,xs = data.shape
   y,x = np.indices(data.shape)
   x_arr = np.arange(xs)

   #fit the curvature
   model = md.models.Chebyshev1D(5)
   fit_c = md.fitting.LinearLSQFitter()
   m = fit_c(model, x[data>0], y[data>0])


   #determine the extraction window
   xc = int(0.5*xs)
   f_arr = data[:,xc]
   y_arr = y[:,xc]
   lab, nlab = nd.label(f_arr>0)
   ymask = y_arr * (lab==lab.max())
   y1 = ymask[ymask>0].min()+dy
   y2 = ymask[ymask>0].max()-dy
   print y1, y2

   mask = (y < m(x)-m(xc)+y1) + (y > m(x)-m(xc)+y2)
   data[mask] = 0

   #imshow(data, aspect='auto', origin='lower')
   #plot(x_arr, m(x_arr)-m(xc)+y1, color='red')
   #plot(x_arr, m(x_arr)-m(xc)+y2, color='red')
   #show()
   return data




def extract_spectra(data, dy=10):
   """Extract the spectra from the image and correct 
      for any curvature
   """

   ys,xs = data.shape
   y,x = np.indices(data.shape)

   #fit the curvature
   model = md.models.Chebyshev1D(5)
   fit_c = md.fitting.LinearLSQFitter()
   m = fit_c(model, x[data>0], y[data>0])

   sdata = 0.0 * data
   for i in range(xs):
       dm = m(0) -  m(i) 
       if dm > 0: 
          sy1 = 0
          sy2 = int(ys-dm)
          ty1 = ys - (sy2-sy1) 
          ty2 = ys 
       else:
          sy1 = int(-dm)
          sy2 = ys
          ty1 = 0
          ty2 = sy2-sy1

       sdata[ty1:ty2,i] = data[sy1:sy2,i]

   #cut the data
   y,x = np.indices(sdata.shape)
   y = y * (sdata > 0)
   y1=max(0, y[(y>0)].min()-5)
   y2=min(len(data), y.max() + 5)
   sdata = sdata[y1:y2,:]
   
   return sdata

def ncor(x, y):
    """Calculate the normalized correlation of two arrays"""
    d=np.correlate(x,x)*np.correlate(y,y)
    if d<=0: return 0
    return np.correlate(x,y)/d**0.5


def calc_shift(arc, ref, dc):
    """Calculate the shift between two points. This only calculates
       an integer shift.

    """
    xarr = np.arange(len(ref))
    best_n = 0
    dx = np.median(dc)
    n_list = []
    for x in dc:
        sarc = np.interp(xarr, xarr+x, arc)
        n = ncor(sarc, ref)
        n_list.append(n)
        if n > best_n: 
           best_n = n
           dx = x
    return dx
  
    

def straight_image(obj, arc, yc=None, max_shift=20):
    """Step through the arc image and calculate the 
       needed shift in order to remove any tilt
       in the slits.  Apply this shift to the data.
    """
    
    ys,xs = obj.shape
    if yc is None: yc =int(0.6*ys)
   
    sobj = 0.0 * obj
    sobj[yc,:] = obj[yc,:]
    limit = arc[yc,:].sum()
    dx = 0
    for i in range(0, yc)[::-1]+range(yc,ys):
        if arc[i,:].sum() > 0.5*limit:
           #calculate the shift
           dx = 0
           dc = np.arange(dx-max_shift, dx+max_shift+1)
           dx = calc_shift(arc[i,:], arc[yc,:], dc)
           if dx > 0:
              x1 = dx
              x2 = xs-x1
              sobj[i,x1:] = obj[i,0:xs-x1]
           else:
              x1 = abs(dx)
              x2 = xs-abs(dx)
              sobj[i,0:x2] = obj[i,x1:]
    return  sobj
    


def get_order(data, orders, arc, sp, order, dy=20):
   """Get a 1D spectra from an image
   """
   mask = (orders==order)
   #print data.shape
   #print mask.shape
   obj = data*mask
   arc = arc*mask
   y,x = np.indices(data.shape)
   y = y * mask
   y1=max(0, y[(y>0)].min()-5)
   y2=min(len(data), y.max() + 5)
   print y1, y2
   obj = 1.0 * obj[y1:y2,:]


   sobj = get_onlyobject(obj, 2, dy=dy)

   #this is the slower, more accurate method
   #sobj = extract_spectra(obj)
   #sarc = extract_spectra(arc)
   #rsobj = straight_image(sobj, sarc)
   #sobj = get_onlyobject(rsobj, 2, dy=10)
   #np.save('rsobj', rsobj)
   #rsobj=np.load('rsobj.npy')
   #imshow(rsobj, origin='lower', aspect='auto')
   #show()


   flux = sobj.sum(axis=0)[::-1]
   xarr = np.arange(len(flux))
   w = 1e7*sp.get_wavelength(xarr)
   return w[50:-100], flux[50:-100]


def find_order(sp, w):
    """Givena wavelength determine the order for that wavelength
    """
    lam = 1e7*sp.calc_centralwavelength()
    n = round(sp.order * lam / w)
    return n

def get_info(wavelength, obsmode):
   if wavelength > red_edge:
      order_file = 'R_'
      sp =  HRSModel(camera_name='hrdet', order=64)
   else:
      order_file = 'H_'
      sp =  HRSModel(camera_name='hbdet', order=104)

   if obsmode.count('HIGH'):
      wave_file = order_file+'HR_Arc.fits'
      order_file += 'HR_Order.fits'
      dy = 10
   elif obsmode == 'MEDIUM RESOLUTION':
      wave_file = order_file+'MR_Arc.fits'
      order_file += 'MR_Order.fits'
      dy = 10
   elif obsmode == 'LOW RESOLUTION':
      wave_file = order_file+'LR_Arc.fits'
      order_file += 'LR_Order.fits'
      dy = 5
   return sp, order_file, wave_file, dy


