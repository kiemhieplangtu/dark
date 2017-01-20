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
xlambda = 0.5

# for x in pl.frange(0., 1., 0.1):
xx = []
yy = []
for j in range(2449,2459):
	alpha  = float(j)*pi/40000.
	r      = 10.
	theta  = 0.
	for i in range(10090):
		x = r*np.cos(theta)
		y = -r*np.sin(theta)

		if (r<1.):
			phi = (float(j) - 1.)*pi/40000.
			break

		dr     = -0.01*np.cos(alpha)
		dtheta = 0.01*np.sin(alpha)/r
		r      = r + dr
		theta  = theta + dtheta
		dalpha = -dtheta - 0.01*xlambda*np.sin(alpha)/r/(r-xlambda)
		alpha  = alpha + dalpha
		xx.append(x)
		yy.append(y)

plt.plot(xx,yy,'r.', ms=1.)
circle = plt.Circle((0, 0), radius=1., fc='k')
plt.gca().add_patch(circle)
plt.xlim(-10.,10.)
plt.ylim(-10.,10.)
plt.grid()
plt.show()