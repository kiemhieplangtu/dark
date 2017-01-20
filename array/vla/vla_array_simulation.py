## IMPORT LIBRARIES ##
import sys, os
import matplotlib.pyplot as plt
import numpy             as np

## CONSTANTS AND PARAMETERS ##
pi     = 2.*np.arccos(0.)
xlamda = 1.e-6                 # Wavelength: lambda = 1mm = 1e-6 km

ns     = 1                     # Number of sources
dec0   = 45.*pi/180.           # Declination of Source RX J0911, dec ~ 5deg50'54'', RA~9h11m27.7s
deltah = 1.*2.*pi/60./24.      # Step for the loop in Hour-angle, 1 min, 2pi/24/60

nhours = 1                     # Hours of observing
nha    = 60*nhours             # hours, time for observing, if want to see snapshot of UV-Coverage -> set nha = 1
ha0    = -deltah*nha/2.        # Initial Hour-angle of the source, Minus sign because I take the clockwise direction
glati  = 34.0784*pi/180.       # Latitude of the array, VLA at latitude 34.0784 degree


## POSITIONS OF SOURCES ##
 # In this simple simulation, I just use a point source with Unit Brightness, so, ns=1
ux  = []  # Position of sources in Sky-coordinates, (l,m)
uy  = []  # Position of sources in Sky-coordinates, (l,m)
dec = []
ha  = []
bi  = []  # Brightness = 1
for i in range(0,ns):
	ux.append(0.)
	uy.append(0.)
	dec.append(dec0)
	ha.append(ha0)
	bi.append(1.) # Brightness = 1

##  POSITIONS OF ANTENNAS, phi=0 IN THE NORTH DIRECTION ##
nanten = 27
rant   = [0.4364, 1.4337, 2.8747, 4.7095, 6.9065, 9.4434, 12.3027, 15.4706, 18.9357,
          0.4840, 1.5899, 3.1881, 5.2229, 7.6595, 10.4728, 13.6438, 17.157, 21.,
          0.484, 1.5899, 3.1881, 5.2229, 7.6595, 10.4728, 13.6439, 17.1572, 21.]  # distances
fi = +5.          
phiant = [0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,0.+fi,
          120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,120.+fi,
          240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi,240.+fi] # Angles

## COORDINATES of ANTENNAS ##
x = np.zeros(nanten) ##
y = np.zeros(nanten) ## 
for i in range(nanten):
	x[i]  = rant[i]*np.cos(phiant[i]*pi/180.)
	y[i]  = rant[i]*np.sin(phiant[i]*pi/180.)

## COORDINATES of BASELINES in UNITS of WAVELENGTH ##
nbase = nanten*(nanten-1)/2
bx    = np.zeros(nbase)
by    = np.zeros(nbase)
nb    = 0                # Infact, nb = nbaseline = nanten*(nanten-1)/2
for i in range(nanten-1):
	for j in range(i+1,nanten):		
		bx[nb] = (x[j] - x[i])/xlamda
		by[nb] = (y[j] - y[i])/xlamda
		nb = nb + 1


## FULFILL the UV-COVERAGE ##
han = ha0
nd  = nanten*(nanten-1)*nha/2
u   = np.zeros(nd)
v   = np.zeros(nd)
nuv = 0            # Number of points in uv-plane
for iha in range(nha):
	han = han + deltah
	for i in range(nb):
		utemp  = -np.cos(han)*by[i]-np.sin(glati)*np.sin(han)*bx[i]
		vtemp  = ( np.sin(glati)*np.sin(dec0)*np.cos(han)+np.cos(glati)*np.cos(dec0) )*bx[i] - np.sin(han)*np.sin(dec0)*by[i]
		elsrc  = np.sin(glati)*np.sin(dec0)+np.cos(glati)*np.cos(dec0)*np.cos(han) # Elevation of source
		if (elsrc > 0):          # check if the source is above the horizon, if below the horizon -> cannot see
			u[nuv] = utemp
			v[nuv] = vtemp
			nuv    = nuv + 1
		if (elsrc < 0):
			print 'Cannot observe'

## PLOT the UV-COVERAGE ##
 # To see the snapshot of uv-coverage, set nha = 1 #	
plt.plot(u*xlamda,v*xlamda, 'b.', ms=2)
plt.plot(-u*xlamda,-v*xlamda, 'b.',ms=2)
plt.title('UV Coverage, VLA, 27 antennas')
plt.grid()
plt.show()

# DIRTY BEAM in SKY-COORDINATES (l,m) ##
 #Write l,m,b to file, then plot later # 
## FOURIER TRANSFRORM, DIRTY-BEAM = FT[UV-SAMPLING] ##
lmax = 512
mmax = 512
for l in range(lmax):
	xl = -0.5 + 0.001953125*l # 1/512 = 0.001953125
	for m in range(mmax):
		xm = -0.5 + 0.001953125*m
		db = 0.
		for k in range(nuv):
			bu = u[k]*xl + v[k]*xm
			db = db + 2.*np.cos( 2.*pi**2*bu/180./3600. ) ## 1rad = pi/180./3600. arcsec

		# print xl,xm,db

## DIRTY IMAGE in SKY-COORDINATES (l,m) ##
 # Write l,m,w to file, then plot later #

## REAL(Vi) and IMG(Vi) at each point in uv-plane ##
vire  = np.zeros(nuv)
viim  = np.zeros(nuv)
phivi = np.zeros(nuv)
for k in range(nuv):
	vire[k] = 0.
	viim[k] = 0.
	for i in range(ns): # Here, ony 1 source, ns = 1
		bu      = u[k]*ux[i] + v[k]*uy[i]  ## ux,uy source position 
		vire[k] = vire[k] + bi[i]*np.cos(2.*pi**2*bu/180./3600.) # =1, Only a Single source at ux = uy = 0, so Real(Vi) = 1, Img(Vi) = 0
		viim[k] = viim[k] + bi[i]*np.sin(2.*pi**2*bu/180./3600.) # =0

	phivi[k] = np.arctan2( viim[k],vire[k] )

## FOURIER TRANSFRORM, DIRTY-MAP = FT[Vi-MAP] ##
lmax = 512
mmax = 512
vi   = np.zeros(nuv)
for l in range(lmax):
	xl = -0.5 + 0.001953125*l
	for m in range(mmax):
		xm = -0.5 + 0.001953125*m
		wi = 0.
		for k in range(nuv):
			bu    = u[k]*xl + v[k]*xm
			vi[k] = np.sqrt( vire[k]**2 + viim[k]**2 )
			wi    = wi + 2.*vi[k]*np.cos( 2.*pi**2*bu/180./3600. - phivi[k] )

		# print xl,xm,wi