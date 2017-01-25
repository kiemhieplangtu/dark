import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from gauss_fit      import gfit

## Get the index of a given velocity #
 #
 # params list v-axis List of Velocity_axis
 # params float vel Velocity
 #
 # return int idx Index of vel in List of velocity_axis
 # 
 # Author Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Get Vrange Indexes ##
 #
 # params 1-D-array v     VLSR
 # params float     lowv  Lower limit
 # params float     upv   Upper limit
 #
 # return list
 #
 # version 01/2017 
 # author Nguyen Van Hiep ##
def get_vrange_id(v, lowv, upv):
	vmax = get_vel_index(v, lowv)
	vmin = get_vel_index(v, upv)
	return [vmin, vmax]

## Multiple (N) Gaussians + offset. ##
 #
 # params list  v    VLSR
 # params float zr   estimated constant zero offset of the data points.
 # params list  h    the array of N estimated heights of the Gaussians.
 # params list  v0   the array of N estimated centers of the Gaussians.
 # params list  w    the array of N estimated halfwidths of the Gaussians.
 #
 # return 1-D-array  tf  The calculated points.
 #
 # version 01/2017 
 # author Nguyen Van Hiep ##
def gcurv(v, zr, h, v0, w):
	#DETERMINE NR OF GAUSSIANS...
	ng = len(h)

	v  = np.array(v)
	h  = np.array(h)
	v0 = np.array(v0)
	w  = np.array(w)

	tf = 0.*v + zr
	for i in range(ng):
		if (w[i] > 0.):
			tf = tf + h[i]*np.exp(- ( (v-v0[i])/(0.6005612*w[i]))**2) # 0.6005612 - 1/e width

	return tf


## ============= MAIN ================ ##
## Class
fit  = gfit()

src  = '3C286'
print 'Fitting...' + src

data = readsav('../idl/gfit_idl_claire/3C286.sav')
# % RESTORE: Restored variable: VLSR.
# % RESTORE: Restored variable: SPEC1.
# % RESTORE: Restored variable: EM_SPEC.
# % RESTORE: Restored variable: SIGMA.
# % RESTORE: Restored variable: EMSIGMA.
# % RESTORE: Restored variable: PSR1.
# % RESTORE: Restored variable: VLSREM.
# % RESTORE: Restored variable: CONT.
# % RESTORE: Restored variable: RMS.

## ABS
vlsr    = data.vlsr[:,0]
tau     = data.spec1[:,0,0]
sigtau  = data.sigma

## EM - Set the velocity range for the emission spectra (to trim Carl's data)
vlsrem  = data.vlsrem
vmin, \
vmax    = get_vrange_id(vlsrem, -100., 75.)
Te      = data.em_spec

## ABS line
plt.plot(vlsr,tau, 'b-', linewidth=2, label='data, Absorption line')
plt.title(src, fontsize=30)
plt.ylabel(r'$\tau$', fontsize=35)
plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
# plt.xlim(0.0, 2.0)
# plt.xlim(-1.0, 6.0)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper left', fontsize=18)
plt.show()

## EM line
plt.plot(vlsrem,Te, 'b-', linewidth=2, label='data, Emission line')
plt.plot(vlsrem[vmin:vmax],Te[vmin:vmax], 'r-', linewidth=2, label='data, Emission line')
plt.title(src, fontsize=30)
plt.ylabel(r'T_{em} [K]', fontsize=35)
plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
# plt.xlim(0.0, 2.0)
# plt.xlim(-1.0, 6.0)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper left', fontsize=18)
plt.show()

# Retrieve initial Gaussian guesses for the absorption spectrum
zro0   = 0.0
hgt0   = [0.0065, 0.0066, 0.009]
cen0   = [-28.48, -14.4, -7.3]
wid0   = [2.32,    4.66,   5.0]
look   = 0

nrg    = len(hgt0)
zro0yn = 0
hgt0yn = [1]*nrg
cen0yn = [1]*nrg
wid0yn = [1]*nrg
corder = 'no'

## WNM
tspin0 = [0.]*nrg
order0 = list(range(nrg))

zro0yn   = 0
tau0yn   = [1]*nrg
cenc0yn  = [1]*nrg
wid0yn   = [1]*nrg
tspin0yn = [0]*nrg

zrownm = 1.
hgtwnm = [0.003]
cenwnm = [-10.4]
widwnm = [21.46]
fwnm   = [0.5]

zrownmyn = 1
hgtwnmyn = 0
cenwnmyn = 0
widwnmyn = 0
fwnmyn   = 0

## Fit these guesses
# tau0  = fit.gcurv(vlsr, zro0, hgt0, cen0, wid0)
# tfit0 = np.exp(-tau0)

## ABS line - From Guesses
# plt.plot(vlsr, tfit0, 'b-', linewidth=2, label='data, Tau abs line')
# plt.title(src, fontsize=30)
# plt.ylabel(r'$\tau$', fontsize=35)
# plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
# # plt.xlim(0.0, 2.0)
# # plt.xlim(-1.0, 6.0)
# plt.grid(True)
# plt.tick_params(axis='x', labelsize=18)
# plt.tick_params(axis='y', labelsize=15)
# plt.legend(loc='upper left', fontsize=18)
# plt.show()

# tfita, sigma, \
# zro1, hgt1, cen1, wid1,\
# sigzro1, sighgt1, sigcen1, sigwid1,\
# cov, problem,\
# nparams = fit.abfit(look, vlsr, tau, [0, len(tau)-1], zro0, hgt0, cen0, wid0, zro0yn, hgt0yn, cen0yn, wid0yn)

tfita, sigma, \
zro1, hgt1, cen1, wid1, tspin1, \
sigzro1, sighgt1, sigcen1, sigwid1, sigtspin1, \
zrownm1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
sigzrownm1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
cov, problem, nloop, \
tb_cont, tb_wnm_tot, tb_cnm_tot, \
exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, vlsr, tau, [0, len(tau)-1],\
	zro0, hgt0, cen0, wid0, tspin0, order0,\
	zro0yn, hgt0yn, cen0yn, wid0yn, tspin0yn, \
	zrownm, hgtwnm, cenwnm, widwnm, fwnm, \
	zrownmyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn, halfasseduse=0.5)

print 'Absorption line: problem...', problem

## ABS line fit & residuals
fig1 = plt.figure(1)
frame1=fig1.add_axes((.1,.3,.8,.7))
plt.plot(vlsr, tfita, 'r-', linewidth=2, label='data, fit')
plt.plot(vlsr,tau, 'b-', linewidth=2, label='data, Absorption line')
plt.ylabel(r'$\tau$', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.tick_params(axis='x', labelsize=18, labelbottom='off')
plt.legend(loc='upper left', fontsize=18)
plt.grid()

frame2=fig1.add_axes((.1,.1,.8,.2))
difference = tau - tfita
plt.plot(vlsr, difference, 'r-', linewidth=2, label='')
frame2.set_ylabel('$Residual$',fontsize=20)
plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
# plt.legend(loc='upper left', fontsize=18)
plt.axhline(y=sigma, xmin=-100, xmax=100, c='k', ls='-.', linewidth=3)
plt.axhline(y=-sigma, xmin=-100, xmax=100, c='k', ls='-.', linewidth=3)
plt.grid()
plt.show()

print '1. sigma ', sigma
print '2. Zro ', zro1
print '3. tau ', hgt1
print '4. cen ', cen1
print '5. wid ', wid1

print ''
print ''
# zrownm1   = 0.0
# hgtwnm1   = [1,1]
# cenwnm1   = [-5,-20]
# widwnm1   = [10,10]
# look      = 0
# nrgwnm    = n_elements(hgtwnm1)
# zrownm1yn = 1
# hgtwnm1yn = 1+intarr(nrgwnm)
# cenwnm1yn = 1+intarr(nrgwnm)
# widwnm1yn = 1+intarr(nrgwnm)
# fwnm1     = fwnm1
# fwnm1yn   = intarr(nrgwnm)