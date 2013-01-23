import sys, os
import numpy as np

from matchstars import linearfunc

import pylab as pl


if __name__=='__main__':

   reffile=sys.argv[1]
   mcol=7
   xa, ya, ma, sa = np.loadtxt(reffile, usecols=(1,2,mcol, 17), unpack=True)

   infile=open('list.transform').readlines()
   tout=open('list.normflux', 'w')

   nobj=len(xa)
   nframe=len(infile)
   m_arr=np.zeros((nobj, nframe))
   s_arr=np.zeros((nobj, nframe))

   m_arr[m_arr==0]=np.nan
   s_arr[m_arr==0]=np.nan

   img_list=[]
   for i, info in enumerate(infile):
       l=info.split()
       img=l[0].strip()
       img_list.append(img)
       t_x=[float(x) for x in l[1:4]]
       t_y=[float(x) for x in l[4:7]]

       txtfile=img.replace('fits', 'cat')
       print txtfile
       xb, yb, mb, sb = np.loadtxt(txtfile, usecols=(1,2,mcol, 17), unpack=True) 
       print xb, yb
       xw=linearfunc(t_x, xb, yb)
       yw=linearfunc(t_y, xb, yb)

       print img
       dlimit=5
       for j in range(len(xa)):
           d=(xw-xa[j])**2+(yw-ya[j])**2
           k=d.argsort()[0]
           if d[k] < dlimit:
              print xa[j], xw[k], ya[j], yw[k], d[k]
              m_arr[j, i]=mb[k]
              s_arr[j, i]=sb[k]

   pl.figure()
   fa=10**(-0.4*ma)
   norm=np.ones(nframe)
   for i in range(nframe):
       fb=10**(-0.4*m_arr[:,i])
       mask=(fb>0)
       n=(fb/fa)[mask].mean()
       norm[i]=-2.5*np.log10(n)
       tout.write('%s  %f %f\n' % (img_list[i], (fb/fa)[mask].mean(),  np.median(s_arr[:,i][mask])/2.354))
   tout.close()


   star_mag=np.zeros(nobj)
   star_std=np.zeros(nobj)

   for j in range(len(xa)):
       cor=m_arr[j,:]-norm
       #pl.plot(cor, marker='o')
       #pl.plot(m_arr[j,:], marker='o')
       star_mag[j]=cor[cor>0].mean()
       star_std[j]=cor[cor>0].std()
   pl.plot(star_mag, star_std, ls='', marker='o')
   pl.show()
