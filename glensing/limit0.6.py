import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import pylab             as pl
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

pi      = 2.*np.arccos(0.)
xlambda = 0.6

# xlamda, alpha_limit
# 0.6 0.219911485751
alim    = 0.219911485751
astep   = 40000
nav     = alim*astep/pi
nav     = int(round(nav,0))
n1      = nav + 100
n2      = nav + 120

## Step of light-ray, ds = 0.01, 1% of Lens's radius
ds = 0.01

# for x in pl.frange(0., 1., 0.1):
xx = []
yy = []
for j in range(n1,n2):
	alpha  = float(j)*pi/astep
	r      = 10.
	theta  = 0.
	for i in range(2500):
		x = r*np.cos(theta)
		y = r*np.sin(theta)

		if (r<1.):
			phi = (float(j) - 1.)*pi/200.
			break

		dr     = -ds*np.cos(alpha)
		dtheta = -ds*np.sin(alpha)/r
		r      = r + dr
		theta  = theta + dtheta
		dalpha = -dtheta - ds*xlambda*np.sin(alpha)/r/(r-xlambda)
		alpha  = alpha + dalpha
		xx.append(x)
		yy.append(y)

plt.plot(xx,yy,'b.', ms=1.)
circle = plt.Circle((0, 0), radius=1., fc='k')
plt.gca().add_patch(circle)
plt.xlim(-4.,10.)
plt.ylim(-10.,2.)
plt.grid()
plt.show()