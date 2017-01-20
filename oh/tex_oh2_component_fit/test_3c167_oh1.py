import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot
from gauss_fit_3c18 import gfit

## Create a line ##
 #
 # params list x x-data
 # params list y y-data
 # params string label Label of line
 # params dict prop Properties of line
 # return dict ret All infor of line
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def cal_bg(vlsr, stock):
	s     = 0.
	count = 0
	for i in range(0,2048):
		if (((vlsr[i] > -35.) and (vlsr[i] < -20.)) or ((vlsr[i] > 20.) and (vlsr[i] < 35.))):
			s     = s+stock[i]
			count = count + 1

	return s/count

## Create a line ##
 #
 # params list x x-data
 # params list y y-data
 # params string label Label of line
 # params dict prop Properties of line
 # return dict ret All infor of line
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def correct_ctrl_chnl(tb):
	for i in range(0,101):
		tb[i][1023] = (tb[i][1022] + tb[i][1024])/2.

	return tb

## Get the index of a given velocity ##
 #
 # params list v_axis Velocity axis
 # params float vel Value of velocity
 #
 # return int idx Index of Velocity in the Vel_list
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Create multiple Gaussian functions ##
 #
 # params list x x axis
 # params list params Parameters
 #
 # return array y Multiple-Gaussian functions
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amp = params[i]    	
        ctr = params[i+1]
        wid = params[i+2]
        y   = y - amp * np.exp( -((x - ctr)/wid)**2)
    return 1.-y


## Read peaks Info ##
 #
 # params string fname File-name
 #
 # return dict ret Info of all peaks
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def peak_info(fname=''):
	ret = {}

	cols = ['idx','src','amp','tau0','wid']
	fmt  = ['i','s','f','f','f']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read()

	for i in range(0,len(info['src'])):
		if info['src'][i] not in ret.keys():
			ret[info['src'][i]] = {}
			k = 0
			ret[info['src'][i]][str(k)] = [info['amp'][i],info['tau0'][i],info['wid'][i]]
		else:
			k = k+1
			ret[info['src'][i]][str(k)] = [info['amp'][i],info['tau0'][i],info['wid'][i]]

	return ret

## Compute the velocity-ranges within FWHM ##
 #
 # params list popt     Fit-result parameters 
 #
 # return List intvl    Velocity-ranges within 2sigma
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_vel_range_fwhm(popt):
	intvl = []
	for j in range(0, len(popt), 3):
		ctrl_vel = popt[1+j]                                # Line center
		sigma    = popt[2+j]/np.sqrt(2)                     # Sigma of gaussian line
		fwhm     = 2.35482*sigma                            # FWHM
		intvl    += [ctrl_vel-fwhm/2., ctrl_vel+fwhm/2.]    # vertical lines at 2*sigma

	return intvl

## retreive a SINGLE value of 408 t_b from haslam et al. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_tb_408(ell,bee,tb_408):
	iell= round(( modangle(ell)/360.)* 1080.)
	ibee= round( (( modangle( bee, 360., negpos=True)+90.)/180.)* 540)

	return tb_408[ibee, iell]

## Convert angle to a specified range by finding the angle modulo the extent of the range. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def modangle(angle, extent=360., negpos=False):
	offset = 0.
	if(negpos):
		offset = extent/2.

	return ((((angle-offset) % extent) + extent) % extent) - offset	

## Assign components of dict ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def dassign(d,wh,val,n=0):
	if(type(val) is list):
		value = val
	else:
		value = []
		for i in range(n):
			value.append(val)

	d[wh] = value

	return d

## Calculate multiple (N) Gaussians + offset. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def gcurv(xdata, zro1, hgt1, cen1, widfit):
	#DETERMINE NR OF GAUSSIANS...
	ngaussians = 1 # cho vui

	tfit = 0.*xdata + zro1
	if (widfit > 0.):
		tfit = tfit + hgt1*np.exp(- ( (xdata-cen1)/(0.6005612*widfit))**2)

	return tfit

## Read vel-range to calculate background ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_bins_to_cal_bg(fname='../sub_data/bins_to_cal_bg.txt'):
	cols     = ['idx','src','avmin1','avmax1','avmin2','avmax2','evmin1','evmax1','evmin2','evmax2']
	fmt      = ['i','s','f','f','f','f','f','f','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()

	return vel_info

## Read vel-range to calculate background ##
 #
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def baseline_from_linear_fit(xd, td, vmin1, vmax1, vmin2, vmax2):
	tb = []
	v  = []
	count = 0
	for i in range(0,2048):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			tb.append(td[i])
			v.append(xd[i])

 	slope,intercept = np.polyfit(v, tb, 1)

 	return intercept

## Compute Tex for 1665 line ##
 #
 # params dict data Data
 # params dict inf408  Info about the Tb_background at 408MHz
 #
 # return fit-results
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tex(data, inf408):
	fit  = gfit()
 	src_list = list(data.la.srcname)
 	src = '3C167'
	n    = src_list.index(src)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	oh_f1  = data.la.cfr_bd1
	vlsr1  = data.la.vlsr_bd1
	oh_f2  = data.la.cfr_bd2
	vlsr2  = data.la.vlsr_bd2

	em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	# Vrange infor #
	fname    = '../data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']	

	## Background ##
	tbg1665 = 2.8+get_tb_408(ell[n], bee[n], inf408.tb_408)*(408./1666.)**2.8 # Tbg from 408MHz
	
	xmin = vmin[n]
	xmax = vmax[n]

	tbg = tbg1665

	# VLSR #
	x = vlsr2[n]
	velplotmin = xmin
	velplotmax = xmax
	# Get Index of the Velocity-range #
	xmax_id  = get_vel_index(x, xmin)   # xmax_index
	xmin_id  = get_vel_index(x, xmax)   # xmin_index
	num_chnl = xmax_id-xmin_id          # Total number of bins
	vrange   = [xmin_id, xmax_id]
	vrange   = [730,1023]

	# Linear fit the baselines of spectra #
	# Vrange infor to calculate baseline - Linear fit #
	bins   = read_bins_to_cal_bg()
	avmin1 = bins['avmin1']
	avmax1 = bins['avmax1']
	avmin2 = bins['avmin2']
	avmax2 = bins['avmax2']

	evmin1 = bins['evmin1']
	evmax1 = bins['evmax1']
	evmin2 = bins['evmin2']
	evmax2 = bins['evmax2']

	tc1665     = baseline_from_linear_fit(x, ab_avg2[n], avmin1[n], avmax1[n], avmin2[n], avmax2[n])
	bg_off1665 = baseline_from_linear_fit(x, em_avg2[n], evmin1[n], evmax1[n], evmin2[n], evmax2[n])
	print '***********************'
	print 'Tc: ',tc1665
	print 'Emission baseline: ',bg_off1665
	print 'Tbg_from_408: ',tbg1665

	print n, '   ', data.la[n].srcname, data.la[n].ell, data.la[n].bee

	## Fit or not #
	zrolvyn = 1
	tauyn   = [1]
	v0yn    = [1]
	widyn   = [1]

	cont  = tc1665
	zrolv = 0.
	tau   = [0.002]
	v0    = [18.2]
	wid   = [2.]
	tex   = [0.]

	zrolvyn = 0
	tauyn   = [1]
	v0yn    = [1]
	widyn   = [1]
	texyn   = [0]
	contyn  = 1

	##  1665  bd1 ###
	xd = vlsr2[n]
	td = ab_avg2[n]

	xd_bd1 = xd
	td_bd1 = td

	tb_tot = fit.tb_exp(xd,zrolv, tau, v0, wid, tex, cont)

	plt.plot(xd, td)
	plt.plot(xd, tb_tot, 'g-')
	plt.title('1665 Initial Fit')
	plt.xlim(xmin, xmax)
	plt.show()
	
	# Fit Absorption line
	tfita, sigma, \
	zrolvfit1, taufit1, v0fit1, widfit1, texfit1, \
	zrolv_er1, tau_er1, v0_er1, wid_er1, tex_er1, \
	contfit1,\
	cont_er1,\
	cov1, nloop1, nloopmax = fit.fit( xd, td, vrange, \
		zrolv, tau, v0, wid, tex, \
    	zrolvyn, tauyn, v0yn, widyn, texyn, \
    	cont,contyn)

	tb_tot = fit.tb_exp(xd, \
		zrolvfit1, taufit1, v0fit1, widfit1, texfit1,\
		contfit1)

	plt.plot(xd,td)
	plt.plot(xd, tb_tot, 'r-')
	plt.title('1665 Absortion line Fit')
	plt.xlim(xmin, xmax)
	plt.show()
	# ===============================#

	# ======== Emission Line =================
	tde_bd1 = em_avg2[n]
	xd      = xd_bd1
	tde     = tde_bd1

	tex   = [3.]
	texyn = [1]

	cont_em    = tbg
	tbaseline  = round(bg_off1665,3)
	tdfit      = tde - tbaseline + tbg

	zrolv = zrolvfit1
	tau   = taufit1
	v0    = v0fit1
	wid   = widfit1
	nrcnm = len(taufit1)

	zrolvyn  = 0
	contemyn = 0
	tauyn    = [0]*nrcnm
	v0yn     = [0]*nrcnm
	widyn    = [0]*nrcnm

	tb_tot = fit.tb_exp(xd, \
			zrolv, tau, v0, wid, tex, cont_em)

	plt.plot(xd,tdfit)
	
	plt.plot(xd,tb_tot,'g-')
	plt.title('1665 Emission First glance')
	plt.xlim(xmin, xmax)
	plt.show()

	tfita, sigma, \
	zrolvfite1, taufite1, v0fite1, widfite1, texfite1, \
	zrolv_ere1, tau_ere1, v0_ere1, wid_ere1, tex_ere1, \
	cont_eme1,\
	cont_ere1,\
	cove1, nloope1, nloopmax = fit.fit( xd, tdfit, vrange, \
		zrolv, tau, v0, wid, tex, \
    	zrolvyn, tauyn, v0yn, widyn, texyn, \
    	cont_em,\
    	contemyn)


	plt.plot(xd,tdfit)
	plt.plot(xd, tb_tot,'r-')
	plt.title('1665 Emission Demo fit')
	plt.xlim(xmin, xmax)
	plt.show()

	print '1665 Ex Demo temps:'
	print texfite1
	print tex_ere1
	print 'baseline linear fit:', tbaseline
	print nloope1

	# sys.exit()

	# bsl = []
	# er  = []
	# for x in np.arange(74.56,74.62,0.001):
	# 	tfita, sigma, \
	# 	zrolvfite1, taufite1, v0fite1, widfite1, texfite1, \
	# 	zrolv_ere1, tau_ere1, v0_ere1, wid_ere1, tex_ere1, \
	# 	cont_eme1,\
	# 	cont_ere1,\
	# 	cove1, nloope1, nloopmax = fit.fit( xd, tde - x + tbg, vrange, \
	# 		zrolv, tau, v0, wid, tex, \
	#     	zrolvyn, tauyn, v0yn, widyn, texyn, \
	#     	cont_em,\
	#     	contemyn)
	# 	bsl.append(x)
	# 	er.append(sigma)

	# tb_tot = fit.tb_exp(xd, \
	# 		zrolvfite1, taufite1, v0fite1, widfite1, texfite1, cont_eme1)

	# plt.plot(xd,tdfit)
	# plt.plot(xd, tb_tot,'r-')
	# plt.title('1665 Emission fit')
	# plt.xlim(xmin, xmax)
	# plt.show()

	# plt.plot(bsl,er)
	# plt.title('Error vs baseline')
	# plt.show()

	## After check Chi2 vs baseline ##
	tbaseline  = 74.571
	tdfit      = tde - tbaseline + tbg

	tfita, sigma, \
	zrolvfite1, taufite1, v0fite1, widfite1, texfite1, \
	zrolv_ere1, tau_ere1, v0_ere1, wid_ere1, tex_ere1, \
	cont_eme1,\
	cont_ere1,\
	cove1, nloope1, nloopmax = fit.fit( xd, tdfit, vrange, \
		zrolv, tau, v0, wid, tex, \
    	zrolvyn, tauyn, v0yn, widyn, texyn, \
    	cont_em,\
    	contemyn)


	plt.plot(xd,tdfit)
	plt.plot(xd, tb_tot,'r-')
	plt.title('1665 Emission Final fit')
	plt.xlim(xmin, xmax)
	plt.show()

	print '1665 Ex Final temps:'
	print texfite1
	print tex_ere1
	print 'baseline final fit:', tbaseline
	print nloope1

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408
cal_tex(data, inf408)

sys.exit()