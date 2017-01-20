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
dec0   = 45.*pi/180.            # Declination of Source RX J0911, dec ~ 5deg50'54'', RA~9h11m27.7s
deltah = 1.*2.*pi/60./24.     # Step of 20 min, 2pi/24/60

nha    = 60*1                     ## hours, time for observing
ha0    = -deltah*nha/2.        ## Hour angle of the source
glati  = 34.0784*pi/180.           # Latitude of observer
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
rant   = [0.4364, 1.4337, 2.8747, 4.7095, 6.9065, 9.4434, 12.3027, 15.4706, 18.9357,
          0.4840, 1.5899, 3.1881, 5.2229, 7.6595, 10.4728, 13.6438, 17.157, 21.,
          0.484, 1.5899, 3.1881, 5.2229, 7.6595, 10.4728, 13.6439, 17.1572, 21.]  # distances
fi = +5.          
phiant = [0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,
          120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,
          240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi] # Angles

nanten = 27
x = np.zeros(nanten)
y = np.zeros(nanten)
for i in range(nanten):
	xtemp = rant[i]*np.cos(phiant[i]*pi/180.)
	ytemp = rant[i]*np.sin(phiant[i]*pi/180.)
	x[i]  = xtemp
	y[i]  = ytemp


nb = 0
npair = nanten*(nanten-1)/2
bx = np.zeros(351)
by = np.zeros(351)
for i in range(nanten-1):
	for j in range(i+1,nanten):		
		bx[nb] = (x[j] - x[i])/xlamda
		by[nb] = (y[j] - y[i])/xlamda
		nb = nb + 1

# for i in range(len(bx)):
# 	print i,bx[i],by[i]

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
		elvqso = np.sin(glati)*np.sin(dec0)+np.cos(glati)*np.cos(dec0)*np.cos(han)
		if (elvqso > 0):
			u[nuv] = utemp
			v[nuv] = vtemp
			nuv = nuv + 1
		if (elvqso < 0):
			print 'aaaa'
		
# plt.plot(u*xlamda,v*xlamda, 'b.', ms=2)
# plt.plot(-u*xlamda,-v*xlamda, 'b.',ms=2)

# # plt.axis('off')
# # plt.savefig('pict.png', bbox_inches='tight', pad_inches = 0)

# # plt.grid()
# plt.show()

# phib = np.zeros(nuv)
# for k in range(nuv):
# 	phib[k] = np.arctan2( u[k],v[k] )

lmax = 512
mmax = 512
for l in range(lmax):
	xl = -0.5 + 0.001953125*l
	for m in range(mmax):
		xm = -0.5 + 0.001953125*m
		db = 0.
		for k in range(nuv):
			bu = u[k]*xl + v[k]*xm
			db = db + 2.*np.cos( 2.*pi**2*bu/180./3600. )

		print xl,xm,db

# x = np.asarray(x, dtype=np.float32)
# y = np.asarray(y, dtype=np.float32)
# plt.plot(x,y, 'r*', ms=5)
# # plt.xlim(-3.,3.)
# # plt.ylim(-3.,3.)
# plt.grid()
# plt.show()

plt.plot(u*xlamda,v*xlamda, 'b.', ms=2)
plt.plot(-u*xlamda,-v*xlamda, 'b.',ms=2)

# plt.axis('off')
# plt.savefig('pict.png', bbox_inches='tight', pad_inches = 0)

plt.grid()
plt.show()

