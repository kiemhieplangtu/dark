import sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import pylab             as pl
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

pi       = 2.*np.arccos(0.)
xlambda  = 0.8
xlambdas = [0., 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 0.999]

## Step of light-ray, ds = 0.01, 1% of Lens's radius
ds = 0.01

for k in range(len(xlambdas)):
	xlambda = xlambdas[k]
	xx = []
	yy = []
	for j in range(0,50):
		alpha  = float(j)*pi/200.
		r      = 10.
		theta  = 0.
		for i in range(2000):
			x = r*np.cos(theta)
			y = r*np.sin(theta)

			if (r<1.):
				phi = (float(j) - 1.)*pi/200.
				# print j,i,phi
				break

			dr     = -ds*np.cos(alpha)
			dtheta = -ds*np.sin(alpha)/r
			r      = r + dr
			theta  = theta + dtheta
			dalpha = -dtheta - ds*xlambda*np.sin(alpha)/r/(r-xlambda)
			alpha  = alpha + dalpha
			xx.append(x)
			yy.append(y)

	# print 'xlamda, alpha_limit'
	# print xlambda, phi+pi/200.

	plt.subplot(3, 3, k+1)
	plt.plot(xx,yy,'b.', ms=1., label='$\lambda$ = R*/R = ' + str(xlambda) +'\n R*: Schwarzchild radius, $2MG/c^{2}$' )
	circle = plt.Circle((0, 0), radius=1., fc='k')
	plt.gca().add_patch(circle)
	plt.xlim(-7.,10.)
	plt.ylim(-7.,7.)

	plt.title('lambda = ' + str(xlambda))
	if(k==3): plt.ylabel('y')
	if(k==7): plt.xlabel('x')
	plt.legend(loc='upper right', fontsize=12)
	plt.grid()

plt.show()


# plt.plot(xx,yy,'b.', ms=1.)
# circle = plt.Circle((0, 0), radius=1., fc='k')
# plt.gca().add_patch(circle)
# plt.xlim(-7.,10.)
# plt.ylim(-7.,7.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('lambda = R*/R = ' + str(xlambda))
# plt.grid()
# plt.show()