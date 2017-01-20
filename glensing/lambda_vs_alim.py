import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import pylab             as pl
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

pi = 2.*np.arccos(0.)
# xlambda = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9999]

## Step of light-ray, ds = 0.01, 1% of Lens's radius
ds = 0.01

xlam = []
alim = []
for xlambda in pl.frange(0., 1., 0.1):
	for j in range(0,150):
		alpha  = float(j)*pi/400.
		r      = 10.
		theta  = 0.
		for i in range(2000):
			x = r*np.cos(theta)
			y = r*np.sin(theta)

			if (r<1.):
				phi = (float(j) - 1.)*pi/400.
				break

			dr     = -ds*np.cos(alpha)
			dtheta = -ds*np.sin(alpha)/r
			r      = r + dr
			theta  = theta + dtheta
			dalpha = -dtheta - ds*xlambda*np.sin(alpha)/r/(r-xlambda)
			alpha  = alpha + dalpha

	xlam.append(xlambda)
	alim.append(phi+pi/400.)

	print 'xlamda, alpha_limit'
	print xlambda, phi+pi/200.

plt.plot(xlam,alim,'r+', ms=15.)
# circle = plt.Circle((0, 0), radius=1., fc='k')
# plt.gca().add_patch(circle)
plt.xlim(-0.1,1.1)
plt.ylim(0.,0.5)
plt.xlabel('Lambda = R*/R, Schwarzchild Radius')
plt.ylabel('Alpha limit (radian), alpha = initial angle of light-ray')
plt.grid()
plt.show()