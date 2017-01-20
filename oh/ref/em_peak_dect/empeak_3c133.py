import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize  import curve_fit
from scipy.io.idl    import readsav
from numpy           import array
from restore         import restore
from plotting        import cplot
from gauss_fit_3c18  import gfit
from scipy.integrate import quad
from mpfit           import mpfit
import copy

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

## Create tau function ##
 #
 # params list x x axis
 # params list params Parameters
 #
 # return array y Multiple-Gaussian functions
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def tau_func(x, *params):
	y = 0.
	for i in range(0, len(params), 4):
		amp = params[i]
		ctr = params[i+1]
		wid = params[i+2]
		tex = params[i+3]
		y   = y + tex*amp * np.exp( -((x - ctr)/wid)**2)

	return y

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


## Bin data-channel up by N neighbouring  points ##
 #
 # params x x-data
 # params y y-data
 # params n nbins
 #
 # return xdata, ydata
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def bin_up(x,t,nbin=4):
	chnl = 2048/nbin

	xx = np.zeros(chnl, dtype=np.float64)
	yy = np.zeros(chnl, dtype=np.float64)
	indx = 0
	for i in range(0,len(x),nbin):
		xtemp = np.sum(x[i:(i+nbin)])
		ytemp = np.sum(t[i:(i+nbin)])
		xx[indx] = xtemp/nbin
		yy[indx] = ytemp/nbin
		indx =  indx + 1

	return xx,yy

## Baseline fit  ##
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
	for i in range(len(xd)):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			tb.append(td[i])
			v.append(xd[i])

	slope,intercept = np.polyfit(v, tb, 1)
	tde = td-(slope*xd+intercept)

	## ============ ##
	s     = 0.
	count = 0
	tbe   = []
	for i in range(len(tb)):
		tbi = tb[i]-(slope*v[i]+intercept)
		tbe.append(tbi)
		count = count + 1
		s     = s + tbi

	u = s/count

	sigma = 0.
	n     = len(tbe)
	for i in range(n):
		sigma = sigma + ((tbe[i] - u)**2)/(n-1)

 	sigma = np.sqrt(sigma)

 	return xd,tde,slope,intercept,tbe,v,u,sigma

## Model of Absorption line and Emission line ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
# def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
def myfunc(p, fjac=None, x=None, y=None, err=None):
	status  = 0
	ng      = (len(p)-1)/3

	fcn = p[0]+p[1]*x
	for i in range(2, len(p),3):
		fcn = fcn + p[i]*np.exp(- ( (x-p[i+1])/(0.6005612*p[i+2]))**2)

	# plt.plot(x,y)
	# plt.plot(x,fcn)
	# plt.show()
	# sys.exit()

	return [status, (y-fcn)/err]

## Dectect peaks of emission line ##
 #
 # params dict data Data
 # params dict inf408  Info about the Tb_background at 408MHz
 #
 # return fit-results
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def peak_dect(data, inf408):
	fit  = gfit()
 	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	oh_f1  = data.la.cfr_bd1
	vlsr1  = data.la.vlsr_bd1
	oh_f2  = data.la.cfr_bd2
	vlsr2  = data.la.vlsr_bd2

	# Oh1665 #
	em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)
	# Oh1667 #
	em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	# Source #
	src = '3C133'
	n = src_list.index(src)

	vlsr   = vlsr1[n]
	oh_fq  = oh_f1[n]
	em_avg = em_avg1[n]
	ab_avg = em_avg2[n]

	# Vrange infor #
	fname    = '../data/em_vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']

	xmin = vmin[n]
	xmax = vmax[n]

	## Background ##
	tbg1665 = 2.8+get_tb_408(ell[n], bee[n], inf408.tb_408)*(408./1666.)**2.8 # Tbg from 408MHz

	tbg = tbg1665
	continuum_em = tbg

	# VLSR #
	x = vlsr
	# Get Index of the Velocity-range #
	xmax_id  = get_vel_index(x, xmin)   # xmax_index
	xmin_id  = get_vel_index(x, xmax)   # xmin_index
	num_chnl = xmax_id-xmin_id          # Total number of bins
	vrange   = [xmin_id, xmax_id]
	# vrange   = [879,1460]
	dv       = (xmax-xmin)/num_chnl


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

	## Bin up ##
	nbin = 1
	x, em_avg = bin_up(x,em_avg,nbin=nbin)
	xd,tde,slope,intercept,tbe,v,u,sigma = baseline_from_linear_fit(x, em_avg, evmin1[n], evmax1[n], evmin2[n], evmax2[n])

	plt.plot(x,em_avg)
	plt.show()
	
	lguess = [0.,0.0025,\
	0.1,7.1,0.2,\
	0.1,8.1,0.2]
	# guessp  = np.array(lguess, dtype='float64')
	pfix    = [False,False,False,False,False,False,False,False]
	plimd   = [[False,False]]*8
	plims   = [[0.,0.]]*8
	pname   = ['bg','slope','amp1','v01','wid1','amp2','v02','wid2']
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]} 

	parinfo = []
	for i in range(len(lguess)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(lguess)):
		parinfo[i]['value']   = lguess[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]

	## Gaussian fit ##
	em_avg = em_avg-(slope*x+intercept)
	xx = x[xmin_id/nbin:xmax_id/nbin]
	yy = em_avg[xmin_id/nbin:xmax_id/nbin]

	fa = {'x':xx, 'y':yy, 'err':yy*0.001}
	mp = mpfit(myfunc, lguess, parinfo=parinfo, functkw=fa, quiet=False)
	print mp.params
	print mp.perror

	## Plot ##
	p = mp.params
	fit = p[0]+p[1]*xx
	amp1 = p[2]
	amp2 = p[5]
	for i in range(2, len(p), 3):
		fit = fit + p[i]*np.exp(- ( (xx-p[i+1])/(0.6005612*p[i+2]))**2)  

	
	print '**********'
	print 'Peak 1: ',amp1,sigma,amp1>3.*sigma
	print 'Peak 2: ',amp2,sigma,amp2>3.*sigma

	plt.plot(x,em_avg)
	plt.plot(xx,fit,'r-')
	plt.title(src)
	plt.xlim(xmin-2.,xmax+2.)
	plt.ylim(-0.03,0.15)
	# plt.axhline(y=sigma,xmin=-30.,xmax=30., linewidth=1, color='r')
	plt.axhline(y=3.*sigma,xmin=-30.,xmax=30.,linewidth=1, color='g', linestyle='dashed')
	plt.axvline(x=xmin,ymin=0.,ymax=1000., linewidth=1, color='k')
	plt.axvline(x=xmax,ymin=0.,ymax=1000., linewidth=1, color='k')
	plt.grid()
	plt.show()

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408
peak_dect(data, inf408)

sys.exit()