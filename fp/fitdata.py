#
#FITDATA--A class for fitting 1 or more gaussian curves to data
#
#
#
from numpy import *
from fit import *

class fitdata:
   """A class for fitting gaussian curves to a spectral features"""
   
   def __init__(self, wave, data, var=1):
       self.data=data
       self.wave=wave
       self.continuum=median(self.data)
       self.sigcont=1.4826*median(abs(self.data-self.continuum))

       self.data=self.data-self.continuum
       self.var=var

       return

   def runfit(self, func, param, data, x=None):
       return fit(func, param, data, x=x, var=self.var)

   def gauss(self,x):
       return self.height() * exp(-0.5*((x-self.mu())/self.sigma())**2)

   def fitgauss(self, mu, width, max):
       self.mu= Parameter(mu)
       self.sigma=Parameter(width)
       self.height = Parameter(max)
       self.param=[self.mu, self.sigma, self.height]
       self.info=self.runfit(self.gauss, self.param, self.data, x=self.wave)

   def doublegauss(self, x):
       y=self.height() * exp(-0.5*((x-self.mu())/self.sigma())**2)
       y += self.heightb() * exp(-0.5*((x-self.mub())/self.sigmab())**2)
       return y

   def fitdouble(self, mu, width, max, mub, widthb, maxb):
       """Fit a double gaussian.  """

       #set up the parameters
       self.mu= Parameter(mu)
       self.sigma=Parameter(width)
       self.height = Parameter(max)
       self.mub= Parameter(mub)
       self.sigmab=Parameter(widthb)
       self.heightb = Parameter(maxb)
       self.param=[self.mu, self.sigma, self.height, self.mub, self.sigmab, self.heightb]

       #run the fit
       self.info=self.runfit(self.doublegauss, self.param, self.data, x=self.wave)


   def multigauss(self, x):
       y=0
       for i in range(len(self.mulist)):
          y += self.heightlist[i]() * exp(-0.5*((x-self.mulist[i]())/self.sigmalist[i]())**2)
       return y

   def fitmulti(self, n, mu,  width, max):
       self.param=[]
       self.mulist=[]
       self.heightlist=[]
       self.sigmalist=[]

       for i in range(n):
          self.mulist.append(Parameter(mu[i]))
          self.heightlist.append(Parameter(max[i]))
          self.sigmalist.append(Parameter(width[i]))
          self.param.extend([self.mulist[i], self.heightlist[i], self.sigmalist[i]])
       self.info=self.runfit(self.multigauss, self.param, self.data, x=self.wave)

