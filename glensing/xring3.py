import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import pylab             as pl
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

# sigmas is the ratio between the lens radius and the lens-source distance
# sigmao is the ratio between the lens radius and the lens-observer distance
# lambda is the ratio between the Schwarzschild radius of the lens and its radius
# zeta is the angular resolution
# omega is the angle between the lens-observer and the lens-source lines
# all variables above are given in microradians (or ppm)

## Read x,y of ring #
 #
 # params string fname Filename
 #
 # return void
 # 
 # Author Van Hiep ##
def read_ring(fname = 'ring.txt'):
   cols = ['x','y']
   fmt  = ['f','f']
   data = restore(fname, 0, cols, fmt)
   return data.read()


pi = 2.*np.arccos(0.)
sigmas,sigmao,xlambda,zeta = [3.,30.,100.,.2]

xx=[]
yy=[]
omega=-1.
for j in range(1):
   omega=omega+10.
   omega = 0.
   xk=1./(1.+sigmao/sigmas)
   zeta2=zeta**2
   sigmaos=sigmao+sigmas
   strength=2.*xlambda*sigmaos
   ring0=xk*np.sqrt(strength)*sigmao/sigmas

   theta0k=0.5*omega
   dthetak=0.5*np.sqrt(omega**2+4.*strength)
   t1k=(theta0k-dthetak-3.*zeta/xk)
   if(t1k < 0.):
      t1k=0.
   
   t2k=(theta0k+dthetak+3.*zeta/xk)
   if(t2k < sigmas/xk):  ## 1st condition
      continue
   if(t1k < sigmas/xk): ## 1st condition
      t1k=sigmas/xk

   t1ksq=t1k**2
   t2ksq=t2k**2
   dtksq=t2ksq-t1ksq
   n=0.
   for i in range(100000):
      thetak=np.sqrt(t1ksq+np.random.rand()*dtksq)
      phi=2.*pi*np.random.rand()
      cphi=np.cos(phi)
      sphi=np.sin(phi)
      rm=(thetak-strength/thetak)/sigmao
      csi=rm*cphi-omega/sigmao
      eta=rm*sphi
      rtheta=(thetak-t1k)/(t2k-t1k)
      if(csi**2+eta**2 > zeta2):
         # print rtheta
         continue

      ring=xk*thetak*sigmao/sigmas
      xring=ring*cphi
      yring=ring*sphi
      xring1=xring/ring0
      yring1=yring/ring0
      n=n+1.
      xx.append(xring)
      yy.append(yring)

# write(33,*)omega,n
plt.plot(xx,yy,'b.', ms=3.)
#ircle = plt.Circle((0, 0), radius=1., fc='k')
# plt.gca().add_patch(circle)
# plt.xlim(-7.,10.)
# plt.ylim(-7.,7.)
plt.grid()
plt.show()