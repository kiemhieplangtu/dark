import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot
from gfit           import gfit

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

## Read vel-range ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_vrange(fname=''):
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()

	return vel_info

## Read Infor of 26 sources with Tex values ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_26src_info_tex(fname=''):
	cols     = ['src','l','b', 'il', 'ib', 'tbg1665', 'tex1665']
	fmt      = ['s','f','f','f','f','f','f']
	src      = restore(fname, 2, cols, fmt)
	info     = src.read()

	return info

## Compute Tex for 1665 line ##
 #
 # params dict data Data
 # params dict inf408  Info about the Tb_background at 408MHz
 #
 # return fit-results
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def main(data, inf408):
	cfit     = gfit()
 	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee


	oh_f1  = data.la.cfr_bd1
	vlsr1  = data.la.vlsr_bd1

	em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	# Gaussian peaks' info - estimate values of Amp, tau0 and Width #
	peak = peak_info('../data/gauss_1665_peaks.txt')

	# Vrange infor #
	vel_info = read_vrange('../data/vel_range.txt')
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']

	# 26 src with no CO #
	fname    = '../data/26src_no_co.txt'
	cols     = ['idx','src']
	fmt      = ['i','s']

	src_no_co = restore(fname, 2, cols, fmt)
	s26info   = src_no_co.read()
	s26src    = s26info['src']

	for src in src_list:
		if(src != 'T0629+10'):
			continue

		n    = src_list.index(src)		
		xmin = vmin[n]
		xmax = vmax[n]

		if (xmin == 0. and xmax == 0.):
			continue

		# VLSR #
		xd    = vlsr1[n]
		td_ab = ab_avg1[n] # Absorption data of 1665
		td_em = em_avg1[n] # Emission data of 1665

		# Get Index of the Velocity-range #
		xmax_id  = get_vel_index(xd, xmin)   # xmax_index
		xmin_id  = get_vel_index(xd, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id           # Total number of bins
		vrange   = [xmin_id, xmax_id]
		# vrange   = [733,1023]


		## Background ##
		tbg1665    = 2.8+get_tb_408(ell[n], bee[n], inf408.tb_408)*(408./1666.)**2.8 # Tbg from 408MHz
		t_on1665   = td_ab
		t_off1665  = td_em
		tc1665     = cal_bg(xd,t_on1665)
		bg_off1665 = cal_bg(xd,t_off1665)
		trx1665    = bg_off1665 - tbg1665

		#Fitting, guess parameters and Fit multiple Gaussians for e(-tau)#
		guess = []
		for c in peak[src]:
			guess += peak[src][c]

		print n, ' <<<>>>> ', src, ell[n], bee[n]

		cont       = tc1665  # Tc
		contoff    = tbg1665 # quan trong, thay doi la thay doi Tex rat nhieu
		tbaseline  = bg_off1665 # baseline of 1665 Em line
		#This case tbaseline is sensitive to Ts, 74.75 -> change Ts #
		# cont       = 35.20
		# tbaseline  = 71.21

		npeak      = 8

		print 'Tc:', cont
		print 'Radio bg: ', contoff
		print 'Baseline (Offsource): ', tbaseline

		tau   = [0.6/cont,1.2/cont,0.5/cont,4./cont,3./cont,2./cont,1.5/cont,1./cont]
		v0    = [0.19,1.44,1.1,3.6,4.7,6.1,6.9,8.0]
		wid   = [.5,.5,.5,.5,.5,.5,.5,.5]
		tex   = [1.e-6]*npeak

		## Fit or not #
		tauyn  = [1]*npeak
		v0yn   = [1]*npeak
		widyn  = [1]*npeak
		tsyn   = [0]*npeak
		contyn = 1

		tfit, error, \
		taufit, v0fit, widfit, tsfit, \
		tau_er, v0_er, wid_er, ts_er, \
		contfit, cont_er, cov, nloop, \
		nloopmax = cfit.fit(xd,td_ab,vrange,tau,v0,wid,tex,cont,
			tauyn,v0yn,widyn,tsyn,contyn)

		## fit with Tex
		tsyn   = [1]*npeak
		ts     = [3.5, 3.4, 2.2, 4.0, 6.6, 5.3, 3.9, 3.6] # for 2nd Fit 

		cont_em = contoff
		tau    = taufit
		v0     = v0fit
		wid    = widfit
		npeaks = len(taufit)

		contyn = 0
		tauyn  = [0]*npeak
		v0yn   = [0]*npeak
		widyn  = [0]*npeak
		tdat   = td_em - tbaseline + contoff

		ta_fit, error, \
		taufit, v0fit, widfit, tsfit, \
		tau_er, v0_er, wid_er, ts_er, \
		contfit, cont_er, cov, nloop, \
	    nloopmax = cfit.fit(xd,tdat,vrange,tau,v0,wid,ts,cont_em,\
				tauyn,v0yn,widyn,tsyn,contyn) # co the sai thu tu 

		print '1665 spin temps:'
		print tsfit
		print ts_er

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav')    # l_cntr, b_cntr, tb_408
main(data, inf408)

sys.exit()