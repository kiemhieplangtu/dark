import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import pylab             as pl
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

# S: source, L: lens, O: observer
# sigmas is the ratio between the lens radius and the lens-source distance
# sigmao is the ratio between the lens radius and the lens-observer distance
# lambda is the ratio between the Schwarzschild radius of the lens and its radius
# zeta is the angular resolution of detector/observer
# omega is the angle between the lens-observer and the lens-source lines
# theta is the angle between the light-ray and the lens-source line, called Polar angle (like a cone)
# define: k = 1. / (1 + sigmaO/sigmaS)
# all variables above are given in microradians (or ppm)

## Solutions (theta/k):
## (theta/k) = w/2 +/- sqrt(w^2+8*lambda*sigmaOS)/2
## O is position of Detector
## M is the position of bending-light-ray hitting the plane (beta) containing the Detector and Normal to the lens-source line
## rm = [ (theta/k) - 2*lambda*SigmaOS/(theta/k) ]/SigmaO
## In plane (beta), the origine is O, so, coordinates of M are: 
## csi = rm*cos(phi) - w/sigmaO
## eta = rm*sin(phi)
## phi is azimuth angle, angle between light-ray and (SLO) plane


## Define constant & params
pi      = 2.*np.arccos(0.)
sigmas  = 3.   ## ratio between the lens radius and the lens-source distance
sigmao  = 30.  ## ratio between the lens radius and the lens-observer distance
xlambda = 100. ## ratio between the Schwarzschild radius of the lens and its radius
zeta    = 0.2  ## Angular resolution of detector/observer
omega   = 20.  ## angle between the lens-observer and the lens-source lines
omegas  = [0., 5., 7.5, 10., 15., 20., 30., 50., 70.]  ## List: angle between the lens-observer and the lens-source lines

for j in range(len(omegas)):
	omega    = omegas[j]
	xk       = 1./(1.+sigmao/sigmas) ## Define: k = 1. / (1 + sigmaO/sigmaS)
	zeta2    = zeta**2
	sigmaos  = sigmao+sigmas
	strength = 2.*xlambda*sigmaos

	## For perfect alignment, omega = 0., later on, I'll normolize ring-sizes (or Arc-size) by ring0
	ring0   = xk*np.sqrt(strength)*sigmao/sigmas

	## (theta/k) = w/2 +/- sqrt(w^2+8*lambda*sigmaOS)/2
	theta0k = 0.5*omega
	dthetak = 0.5*np.sqrt(omega**2+4.*strength)

	## t1k, t2k are range for (theta/k), let (theta/k) vary in a range and check the conditions, 
	## For a ray to be seen by the observer, two conditions must be satisfied: it must
	## avoid the lens and it must hit O plane within the angular resolution eta of the detector.
	## condition1: theta > SigmaS
	## condition2: csi**2+eta**2 > zeta2
	t1k = (theta0k-dthetak-3.*zeta/xk)
	if(t1k < 0.):
		t1k = 0.

	t2k = (theta0k+dthetak+3.*zeta/xk)
	print "theta0k,dthetak,t1k,t2k,ring0,sigmas/xk"
	print theta0k,dthetak,t1k,t2k,ring0,sigmas/xk

	if(t2k < sigmas/xk): ## 1st condition
		sys.exit()
	if(t1k < sigmas/xk): ## 1st condition
		t1k = sigmas/xk

	t1ksq = t1k**2
	t2ksq = t2k**2
	dtksq = t2ksq-t1ksq

	x = []
	y = []
	for i in range(100000):
		thetak = np.sqrt(t1ksq+np.random.rand()*dtksq)
		phi    = 2.*pi*np.random.rand()
		cphi   = np.cos(phi)
		sphi   = np.sin(phi)
		rm     = (thetak-strength/thetak)/sigmao
		csi    = rm*cphi-omega/sigmao
		eta    = rm*sphi

		if( (csi**2+eta**2) > zeta2): ## 2nd condition
			continue

		ring   = xk*thetak*sigmao/sigmas
		xring  = ring*cphi
		yring  = ring*sphi
		xring1 = xring/ring0 ## Nomarlized to ring0, w=0, perfect alignment
		yring1 = yring/ring0 ## Nomarlized to ring0, w=0, perfect alignment

		x.append(xring1)
		y.append(yring1)


	plt.subplot(3, 3, j+1)
	plt.plot(x,y,'b.', ms=5.)
	plt.title('w= ' + str(omega) + ' ppm')
	if(j==3): plt.ylabel('y')
	if(j==7): plt.xlabel('x')
	plt.text(1., 1., 'Detector is moving ->\n w= ' + str(omega) + ' ppm', color='k', fontsize=10)
	plt.xlim(-1.5,2.0)
	plt.ylim(-1.5,1.5)
	plt.grid()

plt.show()










# plt.subplot(3, 3, 2)
# plt.plot(x,y,'b.', ms=1.)
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()

# plt.subplot(3, 3, 3)
# plt.plot(x,y,'b.', ms=1.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()

# plt.subplot(3, 3, 4)
# plt.plot(x,y,'b.', ms=1.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()

# plt.subplot(3, 3, 5)
# plt.plot(x,y,'b.', ms=1.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()

# plt.subplot(3, 3, 6)
# plt.plot(x,y,'b.', ms=1.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()

# plt.subplot(3, 3, 7)
# plt.plot(x,y,'b.', ms=1.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()

# plt.subplot(3, 3, 8)
# plt.plot(x,y,'b.', ms=1.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()

# plt.subplot(3, 3, 9)
# plt.plot(x,y,'b.', ms=1.)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.xlim(-1.5,3.0)
# plt.ylim(-1.5,1.5)
# plt.grid()