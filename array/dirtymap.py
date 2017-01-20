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
dec0   = 5.833*pi/180.            # Declination of Source RX J0911, dec ~ 5deg50'54'', RA~9h11m27.7s
deltah = 20.*2.*pi/60./24.     # Step of 20 min, 2pi/24/60

nha    = 3*6                    ## hours, time for observing
ha0    = -deltah*nha/2.        ## Hour angle of the source
glati  = 45.*pi/180.           # Latitude of observer
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
nanten = 6

rant   = [27.,4.,24.,68.,29.,46.]
phiant = [255.,75.,75.,75.,180.,180.]
x = np.zeros(10)
y = np.zeros(10)
for i in range(nanten):
	xtemp = rant[i]*0.008*np.cos(phiant[i]*pi/180.)
	ytemp = rant[i]*0.008*np.sin(phiant[i]*pi/180.)
	x[i]  = xtemp
	y[i]  = ytemp

# nx = 2
# ny = 2
# nanten = nx*ny
# r0x = 1.5
# r0y = 1.
# for i in range(nx):
# 	for j in range(ny):
# 		x[i+nx*(j-1)] = r0x*(2*i-nx-1)/2.
# 		y[i+nx*(j-1)] = r0x*(2*j-nx-1)/2.


nb = 0
bx = np.zeros(45)
by = np.zeros(45)
for i in range(nanten-1):
	for j in range(i+1,nanten):		
		bx[nb] = ( (x[j] - x[i])/xlamda )
		by[nb] = ( (y[j] - y[i])/xlamda )
		nb = nb + 1

han = ha0
nuv = 0
nd  = nanten*(nanten-1)*nha/2
u   = np.zeros(nd)
v   = np.zeros(nd)
for iha in range(nha):
	han = han + deltah
	for i in range(nb):
		utemp  = -np.cos(han)*by[i]-np.sin(glati)*np.sin(han)*bx[i]
		vtemp  = ( np.sin(glati)*np.sin(dec0)*np.cos(han)+np.cos(glati)*np.cos(dec0) )*bx[i] - np.sin(han)*np.sin(dec0)*by[i]
		elvqso = np.sin(glati)*np.sin(dec0)+np.cos(glati)*np.cos(dec0)*np.cos(han) ## elevation of object
		if (elvqso > 0):
			u[nuv] = utemp
			v[nuv] = vtemp
			# print iha,10+i,u[nuv]*xlamda,v[nuv]*xlamda
			nuv = nuv + 1


# print nuv, len(u),len(v)

vire = np.zeros(3240)
viim = np.zeros(3240)
phivi = np.zeros(3240)
for k in range(nuv):
	vire[k] = 0.
	viim[k] = 0.
	for i in range(ns):
		bu = u[k]*ux[i] + v[k]*uy[i]  ## ux,uy source position 
		vire[k] = vire[k] + bi[i]*np.cos(2.*pi*bu)
		viim[k] = viim[k] + bi[i]*np.sin(2.*pi*bu)

	phivi[k] = np.arctan2( viim[k],vire[k] )

# lmax = 512
# mmax = 512
# vi   = np.zeros(3240)
# for l in range(lmax):
# 	xl = -20.48 + 0.08*l
# 	for m in range(mmax):
# 		xm = -20.48 + 0.08*m
# 		wi = 0.
# 		for k in range(nuv):
# 			bu = u[k]*xl + v[k]*xm
# 			vi[k] = np.sqrt( vire[k]**2 + viim[k]**2 )
# 			wi    = wi + 2.*vi[k]*np.cos( 2.*pi**2*bu/180./3600. - phivi[k] )

# 		print xl,xm,wi

x = np.asarray(x, dtype=np.float32)
y = np.asarray(y, dtype=np.float32)
plt.plot(x,y, 'r^', ms=8)
plt.grid()
plt.show()

plt.plot(u*xlamda,v*xlamda, 'r.')
plt.plot(-u*xlamda,-v*xlamda, 'r.')
plt.grid()
plt.show()

