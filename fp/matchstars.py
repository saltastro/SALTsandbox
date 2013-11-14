import sys, os
import numpy as np

from astrom import triangle

from pyraf import iraf

from scipy import optimize

#Create a text file listing positions of reference star and positions of target stars

#3.  match each image and produce a list of photometry
def old():
    xdiff=np.zeros(len(matchcoo))
    ydiff=np.zeros(len(matchcoo))
    for i in range(len(matchcoo)):
        cooa=matchcoo[i][0]
        coob=matchcoo[i][1]
        xdiff[i]=cooa[0]-coob[0]
        ydiff[i]=cooa[1]-coob[1]
    xb+=np.median(xdiff)
    yb+=np.median(ydiff)
    print np.median(xdiff), np.median(ydiff)
    #now insert the best magnitudes into the list
    for i in range(nobj):
        x1=xa[inda][i]
        y1=ya[inda][i]
        r=((xb-x1)**2+(yb-y1)**2)**2
        j=r.argsort()[0]
        mag_arr[i, k]=mb[j]

def linearfunc(p, x, y):
    return p[0]+p[1]*x+p[2]*y

def residual(p, x1, y1, w):
    f=linearfunc(p, x1, y1)
    return w-f

if __name__=='__main__':

   reffile=sys.argv[1]
 
   mcol=5
   max_nobj=20
   nlim=4

   xa, ya, ma = np.loadtxt(reffile, usecols=(1,2,mcol), unpack=True)
   inda=ma.argsort()

   refcoo=reffile.replace('cat', 'coo')
   fout=open(refcoo, 'w')
   for i in range(len(xa)): fout.write('%5.2f %5.2f\n' % (xa[i], ya[i]))
   fout.close()
   
   infile=open('list.objects').readlines()
   tout=open('list.transform', 'w')
   for img in infile:
       img=img.strip()
       img='f'+img
       txtfile=img.replace('fits', 'cat')
       xb,yb,mb=np.loadtxt(txtfile, usecols=(1,2,mcol), unpack=True)

       txtcoo=txtfile.replace('cat', 'coo')
       fout=open(txtcoo, 'w')
       for i in range(len(xb)): fout.write('%5.2f %5.2f\n' % (xb[i], yb[i]))
       fout.close()
       
       nobj=min(len(xb), max_nobj)
       
       if os.path.isfile('out.coo'): os.remove('out.coo')
       iraf.images.immatch.xyxymatch(txtcoo, refcoo, 'out.coo', 3)
           
       x1,y1,x2,y2=np.loadtxt('out.coo', usecols=(0,1,2,3), unpack=True)

       t_x=[0.,1.,0.]
       t_x=optimize.leastsq(residual, t_x, args=(x2,y2,x1))[0]
       t_y=[0.,1.,0.]
       t_y=optimize.leastsq(residual, t_x, args=(x2,y2,y1))[0]
       for i in range(len(x1)):
           print x2[i], linearfunc(t_x, x1[i],y1[i]),
           print y2[i], linearfunc(t_y, x1[i],y1[i])
       tout.write('%s %f %f %f %f %f %f\n' % (img, t_x[0], t_x[1], t_x[2], t_y[0], t_y[1], t_y[2]))
      
       
  

