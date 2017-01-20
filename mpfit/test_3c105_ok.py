import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import pylab             as pl
import copy

from scipy.io.idl    import readsav
from numpy           import array
from restore         import restore
from gauss_fit_3c18  import gfit

from mpfit           import mpfit
from scipy.integrate import quad


## Correct the central channel ##
 #
 # params 1D-array tb Temperature data
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
 # params int n Order of the Source
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_bins_to_cal_bg(n,fname='../sub_data/bins_to_cal_bg.txt'):
	cols     = ['idx','src','avmin1','avmax1','avmin2','avmax2','evmin1','evmax1','evmin2','evmax2']
	fmt      = ['i','s','f','f','f','f','f','f','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()

	avmin1 = vel_info['avmin1']
	avmax1 = vel_info['avmax1']
	avmin2 = vel_info['avmin2']
	avmax2 = vel_info['avmax2']

	evmin1 = vel_info['evmin1']
	evmax1 = vel_info['evmax1']
	evmin2 = vel_info['evmin2']
	evmax2 = vel_info['evmax2']

	return avmin1[n],avmax1[n],avmin2[n],avmax2[n],evmin1[n],evmax1[n],evmin2[n],evmax2[n]

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
def get_tb(p, x):
	ng  = (len(p)-2)/4
	fcn = p[0]
	for i in range(2, len(p), 4):
		fcn = fcn + p[i]*np.exp(- ( (x-p[i+1])/(0.6005612*p[i+2]))**2)   

	exp_sumtau = np.exp(-fcn)
	tb_cont    = p[1]*exp_sumtau

	# CAL. THE OPACITY TAU
	lenx   = len(x) # sometimes = 2048
	arrtau = np.zeros((lenx, ng),dtype=np.float64)

	for i in range(ng):
		arrtau[:, i] = gcurv(x,p[0],p[i*4+2],p[i*4+3],p[i*4+4])

	sumtau = arrtau.sum(1) #2048

	tb_peaks = np.zeros(lenx)

	# BRIGHT TEMP OF EACH PEAK:
	for i in range(ng):
		temp        = np.reshape(arrtau[:, 0:i+1], (lenx, i+1))
		sumtau_i    = temp.sum(1)
		exp_tau_i   = np.exp(arrtau[:, i] - sumtau_i)
		tb_peaks    = tb_peaks + p[i*4+5] * (1. - np.exp(-arrtau[:,i]) ) * exp_tau_i ## From Carl's paper

	tb_tot = tb_cont + tb_peaks # 2048

	return tb_tot

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
def myabfunc(p, fjac=None, x=None, y=None, err=None):
	status  = 0
	ng      = (len(p)-2)/4

	fcn = p[0]
	for i in range(2, len(p), 4):
		fcn = fcn + p[i]*np.exp(- ( (x-p[i+1])/(0.6005612*p[i+2]))**2)   

	exp_sumtau = np.exp(-fcn)
	tb_cont    = p[1]*exp_sumtau

	return [status, (y-tb_cont)/err]
		
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
def myemfunc(p, fjac=None, x=None, y=None, err=None):
	status  = 0
	ng      = (len(p)-2)/4

	fcn = p[0]
	for i in range(2, len(p), 4):
		fcn = fcn + p[i]*np.exp(- ( (x-p[i+1])/(0.6005612*p[i+2]))**2)   

	exp_sumtau = np.exp(-fcn)
	tb_cont    = p[1]*exp_sumtau

	# CAL. THE OPACITY TAU
	lenx   = len(x) # sometimes = 2048
	arrtau = np.zeros((lenx, ng),dtype=np.float64)

	for i in range(ng):
		arrtau[:, i] = p[0] + p[i*4+2]*np.exp(- ( (x-p[i*4+3])/(0.6005612*p[i*4+4]))**2)

	sumtau = arrtau.sum(1) #2048

	tb_peaks = np.zeros(lenx)

	# BRIGHT TEMP OF EACH PEAK:
	for i in range(ng):
		temp        = np.reshape(arrtau[:, 0:i+1], (lenx, i+1))
		sumtau_i    = temp.sum(1)
		exp_tau_i   = np.exp(arrtau[:, i] - sumtau_i)
		tb_peaks    = tb_peaks + p[i*4+5] * (1. - np.exp(-arrtau[:,i]) ) * exp_tau_i ## From Carl's paper

	tb_tot = tb_cont + tb_peaks # 2048
	return [status, (y-tb_tot)/err]

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
	tfit = 0.*xdata + zro1
	if (widfit > 0.):
		tfit = tfit + hgt1*np.exp(- ( (xdata-cen1)/(0.6005612*widfit))**2)

	return tfit

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

## 1-sigma Error of tb-data ##
 #
 # params xd x-data
 # params td y-data
 # params vmin1 vel-range
 # params vmax1 vel-range
 # params vmin2 vel-range
 # params vmax2 vel-range
 #
 # return sigma
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_tb_sigma(xd, td, vmin1, vmax1, vmin2, vmax2):
	tb = []
	v  = []
	for i in range(0,2048):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			tb.append(td[i])
			v.append(xd[i])

	slope,intercept = np.polyfit(v, tb, 1)

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

 	ery = np.zeros(td.shape)
 	for i in range(2048):
		ery[i] = sigma

 	return ery

## Compute Tex for 1665 line ##
 #
 # params int n order of Source
 #
 # return float xmin, xmax Vel-range to fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def vel_range(n):
	fname    = '../data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']	
	
	xmin = vmin[n]
	xmax = vmax[n]

	return xmin, xmax

## Get genaral info of the source ##
 #
 # params dict data Data
 # params str src Source-name
 #
 # return float xmin, xmax Vel-range to fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_src_info(data,src):
	src_list = list(data.la.srcname)
	n        = src_list.index(src)

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

	return n,ell[n],bee[n],oh_f1[n],vlsr1[n],oh_f2[n],vlsr2[n],em_avg1[n],ab_avg1[n],em_avg2[n],ab_avg2[n]

## Fit the absorption line  ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def ab_fit(xd,td,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2):
	npar    =  len(lguess)
	guessp  = np.array(lguess, dtype='float64')
	plimd   = [[False,False]]*npar
	plims   = [[0.,0.]]*npar
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
	pname   = ['bg']*10
	pfix    = [False]*npar

	count   = 1
	for i in range(npar):
		if(i==0):
			pfix[i] = True

		if(i==1):
			pname[i] = 'Tcont'

		if( ((i-2)%4 == 0) and (i>1) ):
			pname[i] = 'tau'+str(count)

		if( ((i-3)%4 == 0) and (i>1) ):
			pname[i] = 'v0'+str(count)

		if( ((i-4)%4 == 0) and (i>1) ):
			pname[i] = 'wid'+str(count)

		if( ((i-5)%4 == 0) and (i>1) ):
			pname[i] = 'tex'+str(count)
			pfix[i]  = True
			count    = count + 1

	parinfo = []
	for i in range(len(guessp)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		parinfo[i]['value']   = guessp[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]

	##  1665 ###
	erry = get_tb_sigma(xd, td, evmin1, evmax1, evmin2, evmax2)
	x    = xd[xmin_id:xmax_id]
	y    = td[xmin_id:xmax_id]
	erry = erry[xmin_id:xmax_id]

	x    = x.astype(np.float64)
	y    = y.astype(np.float64)
	er   = erry.astype(np.float64)

	fa = {'x':x, 'y':y, 'err':er}
	mp = mpfit(myabfunc, guessp, parinfo=parinfo, functkw=fa, quiet=True)

	## ********* Results ********* ##
	print '********* Results *********'
	abp   = mp.params
	abper = mp.perror
	for i in range(len(parinfo)):
		print "%s = %f +/- %f" % (parinfo[i]['parname'],abp[i],abper[i])
	## Plot ##
	fit   = abp[0]
	for i in range(2, len(abp), 4):
		fit = fit + abp[i]*np.exp(- ( (x-abp[i+1])/(0.6005612*abp[i+2]))**2)  

	fit = np.exp(-fit)
	fit = abp[1]*fit
	plt.plot(x,y)
	plt.plot(x,fit,'r-')
	plt.show()

	return abp,abper,npar,parbase,pname,parinfo

## Fit the absorption line to get Baseline ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def em_fit_for_baseline(xd,tde,\
	abp,abper,npar,parbase,pname,parinfo,\
	base_range,tbg1665,bg_off1665,nbinup,\
	xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2):
	## Bin up ##
	nbin = nbinup
	# xd, td = bin_up(xd,td,nbin=nbin)

	cont_em    = tbg1665
	tbaseline  = round(bg_off1665,3)

	elguess    = [0.]*npar
	plimd      = [[False,False]]*npar
	plims      = [[0.,0.]]*npar
	pfix       = [False]*npar

	for i in range(npar):
		elguess[i] = abp[i]
		if(i==0):
			pfix[i]    = True
			elguess[i] = 0.
			plimd[i]   = [False,False]

		if(i==1):
			pfix[i]    = True
			elguess[i] = cont_em
			plimd[i]   = [False,False]		

		if( ((i-2)%4 == 0) and (i>1) ): # amp
			pfix[i]  = True
			plimd[i] = [False,False]

		if( ((i-3)%4 == 0) and (i>1) ): # v0
			pfix[i]  = False
			plimd[i] = [True,True]
			plims[i] = [abp[i]-abper[i],abp[i]+abper[i]]

		if( ((i-4)%4 == 0) and (i>1) ): # wid
			pfix[i]  = False
			plimd[i] = [True,True]
			plims[i] = [abp[i]-abper[i],abp[i]+abper[i]]

		if( ((i-5)%4 == 0) and (i>1) ): # tex
			elguess[i] = 5.
			pfix[i]    = False
			plimd[i]   = [False,False]

	guessp   = np.array(elguess, dtype='float64')
	eparinfo = []
	for i in range(len(guessp)):
		eparinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		eparinfo[i]['value']   = guessp[i]
		eparinfo[i]['fixed']   = pfix[i]
		eparinfo[i]['parname'] = pname[i]
		eparinfo[i]['limited'] = plimd[i]
		eparinfo[i]['limits']  = plims[i]

	## Do the fit ##
	bsl = []
	erl = []
	for xbsl in np.arange(base_range[0],base_range[1],0.001):
		### X-data, Y-data and 1-sigma Y-Error ###
		tdef = tde - xbsl + tbg1665
		erry = get_tb_sigma(xd, tdef, evmin1, evmax1, evmin2, evmax2)
		xx   = xd[xmin_id/nbin:xmax_id/nbin]
		yy   = tdef[xmin_id/nbin:xmax_id/nbin]
		erry = erry[xmin_id:xmax_id]
		
		xx = xx.astype(np.float64)
		yy = yy.astype(np.float64)
		er = erry.astype(np.float64)
		fa = {'x':xx, 'y':yy, 'err':er}

		mp = mpfit(myemfunc, guessp, parinfo=eparinfo, functkw=fa, quiet=True)

		## Residuals ##
		p      = mp.params
		fit    = get_tb(p,xx)
		resid  = yy - fit
		resid2 = np.square(resid)
		sigsq  = resid2.sum()
		sigma  = sigsq**0.5

		bsl.append(xbsl)
		erl.append(sigma)

	plt.plot(bsl,erl)
	plt.title('Error vs baseline')
	plt.grid()
	plt.show()
	baseline = bsl[erl.index(min(erl))]

	return guessp,pfix,plimd,eparinfo,baseline

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
	fit        = gfit()
 	src        = '3C105'
 	base_range = [52.22,52.26]
 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src)
 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = read_bins_to_cal_bg(n)

 	## Absorption and Emission data ##
 	xd  = vlsr1
	td  = ab_avg1
	tde = em_avg1

	## Background and Vrange ##
	tbg1665    = 2.8+get_tb_408(ell,bee,inf408.tb_408)*(408./1666.)**2.8 # Tbg from 408MHz
	xmin, xmax = vel_range(n)
	tc1665     = baseline_from_linear_fit(xd, ab_avg1, avmin1, avmax1, avmin2, avmax2)
	bg_off1665 = baseline_from_linear_fit(xd, em_avg1, evmin1, evmax1, evmin2, evmax2)
	
	# Get Index of the Velocity-range #
	xmax_id  = get_vel_index(xd, xmin)   # xmax_index
	xmin_id  = get_vel_index(xd, xmax)   # xmin_index
	num_chnl = xmax_id-xmin_id           # Total number of bins
	vrange   = [xmin_id, xmax_id]
	dv       = (xmax-xmin)/num_chnl
	# vrange   = [879,1460]

	

	## Fit Absorption line to obtain Tau, V0, Width ##
	lguess  = [0.,tc1665,\
	0.025,8.1,0.8,0.,\
	0.01,10.4,0.5,0.\
	]
	abp,abper,npar,parbase,pname,parinfo = ab_fit(xd,td,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2)

	# Fit Emission line to find the best baseline #
	nbinup = 1
	guessp,pfix,plimd,eparinfo,baseline = em_fit_for_baseline(xd,tde,\
																abp,abper,npar,parbase,pname,parinfo,\
																base_range,tbg1665,bg_off1665,nbinup,\
																xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2)

	## Final fit after finding the best Baseline, baseline = 52.237##
	tdf  = tde - baseline + tbg1665
	erry = get_tb_sigma(xd, tdf, evmin1, evmax1, evmin2, evmax2)
	xe   = xd[xmin_id/nbinup:xmax_id/nbinup]
	ye   = tdf[xmin_id/nbinup:xmax_id/nbinup]
	erry = erry[xmin_id:xmax_id]
	
	x  = xe.astype(np.float64)
	y  = ye.astype(np.float64)
	er = erry.astype(np.float64)
	fa = {'x':xe, 'y':ye, 'err':er}
	mp = mpfit(myemfunc, guessp, parinfo=eparinfo, functkw=fa, quiet=True)



	
	## ********* Results ********* ##
	print '***********************'
	print n, ' ', data.la[n].srcname, data.la[n].ell, data.la[n].bee
	print 'Tc: ',                tc1665
	print 'Emission baseline: ', bg_off1665
	print 'Tbg_from_408: ',      tbg1665

	print '********* Absorption Results *********'
	for i in range(len(parinfo)):
		print "%s = %f +/- %f" % (parinfo[i]['parname'],abp[i],abper[i])

	print '********* Emission Results *********'
	emp   = mp.params
	emper = mp.perror
	fit   = get_tb(emp, x)
	for i in range(len(eparinfo)):
		print "%s = %f +/- %f" % (eparinfo[i]['parname'],emp[i],emper[i])

	## Plot ##
	plt.plot(x,y)
	plt.plot(x,fit,'r-')
	plt.title('Final fit of Emission line')
	plt.show()




    ## Calculate N(OH) ##
	popt = []
	for i in range(2, len(emp), 4):
		popt += [emp[i], emp[i+1], emp[i+2],emp[i+3]]

	stau_fit = quad(tau_func, xmin,xmax, args=tuple(popt))
	print '0) Tex[i]*integral(tau)'
	print '    ', stau_fit

	noh_fit = 2.39854792704*stau_fit[0]  # x10^14
	print '1) N(OH):'
	print '    ', noh_fit, 'x10^14'

	print '2) Tbaseline:'
	print '    ', baseline

#============== MAIN ==============#
data   = readsav('../oh/data/makelines.sav') #data.la
inf408 = readsav('../oh/data/tb_408.sav') # l_cntr, b_cntr, tb_408
cal_tex(data, inf408)

sys.exit()