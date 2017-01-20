import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

pi     = 2.*np.arccos(0.)
ns     = 1                     # Number of sources
dec0   = 6.*pi/180.            # Declination of Source RX J0911, dec ~ 5deg50'54'', RA~9h11m27.7s
deltah = 20.*2.*pi/60./24.     # Step of 20 min, 2pi/24/60

nha    = 1                     ## hours, time for observing
ha0    = -deltah*nha/2.        ## Hour angle of the source
ha0    = 0.
glati  = 21.*pi/180.           # Latitude of observer
xlamda = 1.e-6                 # lambda = 1mm = 1e-6 km

## Position of the sources ##
ux  = []
uy  = []
dec = []
ha  = []
bi  = []
for i in range(0,ns):
	ux.append(0.)
	uy.append(0.)
	dec.append(dec0)
	ha.append(ha0)
	bi.append(1.)

##  Position of Antennas ##
rant   = [1.,2.,3.,1.,2.,3.,1.,2.,3.]  # distances
phiant = [30.,30.,30.,150.,150.,150.,270.,270.,270.] # Angles

nanten = 9
x = np.zeros(nanten)
y = np.zeros(nanten)
for i in range(nanten):
	xtemp = rant[i]*np.cos(phiant[i]*pi/180.)
	ytemp = rant[i]*np.sin(phiant[i]*pi/180.)
	x[i]  = xtemp
	y[i]  = ytemp


nb = 0
bx = np.zeros(36)
by = np.zeros(36)
for i in range(nanten-1):
	for j in range(i+1,nanten):		
		bx[nb] = (x[j] - x[i])/xlamda
		by[nb] = (y[j] - y[i])/xlamda
		nb = nb + 1

# for i in range(len(bx)):
# 	print i,bx[i],by[i]

han = ha0
nuv = 0
u   = np.zeros(3240)
v   = np.zeros(3240)
for iha in range(nha):
	han = han + deltah
	for i in range(nb):
		utemp  = -np.cos(han)*by[i]-np.sin(glati)*np.sin(han)*bx[i]
		vtemp  = ( np.sin(glati)*np.sin(dec0)*np.cos(han)+np.cos(glati)*np.cos(dec0) )*bx[i] - np.sin(han)*np.sin(dec0)*by[i]
		elvqso = np.sin(glati)*np.sin(dec0)+np.cos(glati)*np.cos(dec0)*np.cos(han)
		if (elvqso < 0):
			print 'aaaa'
		u[nuv] = utemp
		v[nuv] = vtemp

		nuv = nuv + 1

vire = np.zeros(3240)
viim = np.zeros(3240)
phivi = np.zeros(3240)
for k in range(nuv):
	vire[k] = 0.
	viim[k] = 0.
	for i in range(ns):
		bu = u[k]*ux[i] + v[k]*uy[i]
		vire[k] = vire[k] + bi[i]*np.cos(2.*pi*bu)
		viim[k] = viim[k] + bi[i]*np.sin(2.*pi*bu)

	phivi[k] = np.arctan2( viim[k],vire[k] )

lmax = 512
mmax = 512
vi   = np.zeros(3240)
for l in range(lmax):
	xl = -0.5 + 0.001953125*l
	for m in range(mmax):
		xm = -0.5 + 0.001953125*m
		wi = 0.
		for k in range(nuv):
			bu = u[k]*xl + v[k]*xm
			vi[k] = np.sqrt( vire[k]**2 + viim[k]**2 )
			wi    = wi + 2.*vi[k]*np.cos( 2.*pi**2*bu/180./3600. - phivi[k] )

		print xl,xm,wi

x = np.asarray(x, dtype=np.float32)
y = np.asarray(y, dtype=np.float32)
plt.plot(x,y, 'r.')
plt.xlim(-3.,3.)
plt.ylim(-3.,3.)
plt.grid()
plt.show()

plt.plot(u*xlamda,v*xlamda, 'r.')
plt.plot(-u*xlamda,-v*xlamda, 'r.')
plt.grid()
plt.show()

