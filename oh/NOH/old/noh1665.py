import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import copy

from scipy.io.idl        import readsav
from restore             import restore
from mpfit               import mpfit
from scipy.integrate     import quad


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

## Read peaks Info ##
 #
 # params string fname File-name
 # return list guessp Guess-Parameters 
 #        list base_range Range of Baselines (in this script, don't care about it)
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def peak_info(source,fname='../data/gauss_1665_peaks.txt'):
	ret = {}

	cols = ['idx','src','tau','v0','wid','bmin','bmax']
	fmt  = ['i','s','f','f','f','f','f']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read()
	bmin = info['bmin']
	bmax = info['bmax']
	src  = info['src']
	tau  = info['tau']
	v0   = info['v0']
	wid  = info['wid']

	for i in range(0,len(src)):
		if info['src'][i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['base_range'] = [bmin[i],bmax[i]]
			ret[src[i]]['guessp']     = [tau[i],v0[i],wid[i]]
		else:
			ret[src[i]]['base_range'] = [bmin[i],bmax[i]]
			ret[src[i]]['guessp']     += [tau[i],v0[i],wid[i]]

	return ret[source]['guessp'], ret[source]['base_range']

## Create tau function ##
 #
 # params list x x axis
 # params list params Parameters
 #
 # return array y Multiple-Gaussian function
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
		y   = y + tex*amp * np.exp( -((x - ctr)/(0.6005612*wid))**2)
	
	return y

## Retreive a SINGLE value of 408 t_b from haslam et al. ##
 #
 # params float ell Galactic-longitude
 # params float bee Galactic-latitude
 #
 # return float Tbg_408 Background-Temperature at (l,b)
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

## Read Tbg of 408MHz from Healpy map ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def read_tbg408_healpy(fname='../result/bg408_to_compare.txt'):
	cols = ['src','l','b', 'il', 'ib', 'tbg', 'l-idx','b-idx','tbg1','tbg_hpy']
	fmt  = ['s','f','f','f','f','f','f','f','f','f']
	src  = restore(fname, 2, cols, fmt)
	info = src.read()

	src  = info['src']
	gl   = info['l']
	gb   = info['b']
	il   = info['il']
	ib   = info['ib']
	tbg  = info['tbg']
	lid  = info['l-idx']
	bid  = info['b-idx']
	bg1  = info['tbg1']
	bgh  = info['tbg_hpy']

	ret  = {}
	for i in range(len(src)):
		ret[src[i]] = bgh[i]

	return ret

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
 # params 1-D array xd x-data
 # params 1-D array td T-data
 # params float vmin1 Range of velocity
 # params float vmax1 Range of velocity
 # params float vmin2 Range of velocity
 # params float vmax2 Range of velocity
 #
 # return slope and intercept of the Linear fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def baseline_from_linear_fit(xd,td,vmin1,vmax1,vmin2,vmax2,fit=True):
	if(fit):
		tb = []
		v  = []
		count = 0
		for i in range(0,len(xd)):
			if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
				tb.append(td[i])
				v.append(xd[i])
				count = count + 1

	 	slope,bsline = np.polyfit(v, tb, 1)
	else:
	 	bsline = 0.
		count  = 0.
		slope  = 0.
		for i in range(0,len(xd)):
			if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
				bsline = bsline + td[i]
				count  = count + 1

		bsline = bsline/count

 	sigma = 0.
 	for i in range(0,len(xd)):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			sigma = sigma + (td[i]-(slope*xd[i] + bsline))**2

	sigma = np.sqrt(sigma/(count-1))
 	return bsline, sigma

## Model of Absorption line and Emission line ##
 #
 # params 1-D array p Parameters
 # params fjac
 # params 1-D array x x-data
 # params 1-D array y y-data
 # params 1-D arrya err 1-sigma error of y-data
 #
 # return status of the fit
 #        weighted residuals
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
# def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
def myabfunc(p, fjac=None, x=None, y=None, err=None):
	status  = 0
	ng      = (len(p)-1)/3

	fcn = 0.
	for i in range(1, len(p), 3):
		fcn = fcn + p[i]*np.exp(- ( (x-p[i+1])/(0.6005612*p[i+2]))**2)   

	exp_sumtau = np.exp(-fcn)
	tb_cont    = p[0]*exp_sumtau # P[0] is Tc

	return [status, (y-tb_cont)/err]

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

	xx   = np.zeros(chnl, dtype=np.float64)
	yy   = np.zeros(chnl, dtype=np.float64)
	indx = 0
	for i in range(0,len(x),nbin):
		xtemp    = np.sum(x[i:(i+nbin)])
		ytemp    = np.sum(t[i:(i+nbin)])
		xx[indx] = xtemp/nbin
		yy[indx] = ytemp/nbin
		indx     =  indx + 1

	return xx,yy

## 1-sigma Error of tb-data ##
 #
 # params 1-D array xd x-data
 # params 1-D array td y-data
 # params float vmin1 vel-range
 # params float vmax1 vel-range
 # params float vmin2 vel-range
 # params float vmax2 vel-range
 #
 # return sigma
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_tb_sigma(xd, td, vmin1, vmax1, vmin2, vmax2):
	tb = []
	v  = []
	for i in range(len(xd)):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			tb.append(td[i])
			v.append(xd[i])

	v  = np.asarray(v, dtype=np.float64)
	tb = np.asarray(tb, dtype=np.float64)

 	return np.std(tb)

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

# Read Tbg of 1666 from 408MHz #
#
# params string fname Filename
#
# return dict info of Tbg
# 
# Author Van Hiep
##
def read_tex_carl(fname = '../sub_data/tex_oh1165.txt'):

	cols = ['idx','src','amp','v0','wid','ts1','er1','ts2','er2','tbg', 'ts_carl', 'tau']
	fmt  = ['i','s','f','f','f','f','f','f','f','f','f','f']
	data = restore(fname, 3, cols, fmt)
	dat  = data.read()
	src  = dat['src']
	tau  = dat['tau']
	ts   = dat['ts_carl']

	ret  = {}
	for i in range(0,len(src)):
		if src[i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['tau']     = [tau[i]]
			ret[src[i]]['ts_carl'] = [ts[i]]
		else:
			ret[src[i]]['tau']     += [tau[i]]
			ret[src[i]]['ts_carl'] += [ts[i]]

	return ret

## Get genaral info of the source ##
 #
 # params dict data Data
 # params str src Source-name
 #
 # return general Infos of source
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_src_info(data, src, src_list):
	n       = src_list.index(src)

	ra50    = data.la.ra1950
	dec50   = data.la.dec1950
	ell     = data.la.ell
	bee     = data.la.bee

	oh_f1   = data.la.cfr_bd1
	vlsr1   = data.la.vlsr_bd1
	oh_f2   = data.la.cfr_bd2
	vlsr2   = data.la.vlsr_bd2

	em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	return n,ell[n],bee[n],oh_f1[n],vlsr1[n],oh_f2[n],vlsr2[n],em_avg1[n],ab_avg1[n],em_avg2[n],ab_avg2[n]

## Compute the velocity-ranges within FHWM ##
 #
 # params list popt     Fit-result parameters 
 # return List intvl    Velocity-ranges within FHWM
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_peak_vel_range(popt):
	intvl = []
	for j in range(1, len(popt), 3):
		v0    = popt[1+j]                  # Line center
		wid   = popt[2+j]                  # Sigma of gaussian line
		intvl += [v0-wid/2., v0+wid/2.]    # vertical lines at 1/e

	return intvl

## Fit the absorption line  ##
 #
 # params string src Source
 # params 1-D array xd x-data 
 # params 1-D array td y-data 
 # params 1-D array lguess Guess-parameters
 #
 # return Infos and Parameters of the fit
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def ab_fit(src,xd,td,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2):
	npar    =  len(lguess)
	guessp  = np.array(lguess, dtype='float64')
	plimd   = [[False,False]]*npar
	plims   = [[0.,0.]]*npar
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
	pname   = ['bg']*npar
	pfix    = [False]*npar

	if(src == 'T0629+10'):
		pfix[3]  = True
		pfix[6]  = True
		plimd[3] = [True,True]
		plimd[6] = [True,True]
		plims[3] = [guessp[3]-1.,guessp[3]+1.]
		plims[6] = [guessp[6]-1.,guessp[6]+1.]

	count   = 1
	for i in range(npar):
		if(i==0):
			pname[i] = 'Tcont'

		if( ((i-1)%3 == 0) and (i>0) ):
			pname[i] = 'tau'+str(count)

		if( ((i-2)%3 == 0) and (i>0) ):
			pname[i] = 'v0'+str(count)

		if( ((i-3)%3 == 0) and (i>0) ):
			pname[i] = 'wid'+str(count)
			count    = count + 1

	parinfo = []
	for i in range(len(guessp)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		parinfo[i]['value']   = guessp[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]
		parinfo[i]['lims']    = plims[i]

	##  data and 1-sigma uncertainty ##
	x    = xd[xmin_id:xmax_id]
	y    = td[xmin_id:xmax_id]
	erry = np.zeros(y.shape, dtype=np.float64)
	sig  = get_tb_sigma(xd, td, evmin1, evmax1, evmin2, evmax2)
	erry.fill(sig)

	x    = x.astype(np.float64)
	y    = y.astype(np.float64)
	erry = erry.astype(np.float64)

	fa   = {'x':x, 'y':y, 'err':erry}
	mp   = mpfit(myabfunc, guessp, parinfo=parinfo, functkw=fa, quiet=True)

	## Fit values of Tau ##
	abp   = mp.params
	abper = mp.perror
	fit   = 0.
	for i in range(1, len(abp), 3):
		fit = fit + abp[i]*np.exp(- ( (x-abp[i+1])/(0.6005612*abp[i+2]))**2)  

	fit = np.exp(-fit)
	# fit = abp[1]*fit

	return x,fit,y/abp[0],abp,abper,npar,parbase,pname,parinfo

## Compute Tex for 1665 line ##
 #
 # params dict data     Data
 # params dict inf408   Info about the Tb_background at 408MHz
 # params int  bd       OH665 or OH667
 #
 # return Tex and N(OH)
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def cal_tex_print(data,inf408,bd=1):
	bg408    = read_tbg408_healpy()
 	src_list = list(data.la.srcname)

 	for src in src_list:
	 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = read_bins_to_cal_bg(n)
	 	xmin, xmax                                                        = vel_range(n)

	 	# if (src != '4C13.67'):
			# continue
		if (xmin == 0. and xmax == 0.):
			continue

	 	## Absorption and Emission data ##
	 	xd  = vlsr1
		td  = ab_avg1
		tde = em_avg1
		cst = 3.99757843817
		frq = 1665.402
		pfl = '../data/gauss_1665_peaks.txt'
		if(bd == 2):
			xd  = vlsr2
			td  = ab_avg2
			tde = em_avg2
			cst = 2.21841824609
			frq = 1667.359
			pfl = '../data/gauss_1667_peaks.txt'

		# xd,td = bin_up(xd,td,nbin=1)

		## Background ##
		tbg1665    = 2.8+get_tb_408(ell,bee,inf408.tb_408)*(408./frq)**2.8 # Tbg from 408MHz
		tbg1665    = 2.8+bg408[src]*(408./frq)**2.8 # Tbg from 408MHz
		if(src=='3C123'):
			tbg1665 = 26.

		## BACKGROUNDS and THEIR UNCERTAINTIES ##
		tc1665,tc_er        = baseline_from_linear_fit(xd, td, avmin1, avmax1, avmin2, avmax2,fit=False)
		bg_off1665,bgoff_er = baseline_from_linear_fit(xd, tde, evmin1, evmax1, evmin2, evmax2,fit=False)
		trx                 = bg_off1665 - tbg1665

		## 1-SIGMA STANDARD DEVIATION OF Tabspt and Temmission ##
		tab_sigma  = tc_er
		tem_sigma  = bgoff_er
		trx_sigma  = bgoff_er
		ton_sig    = np.sqrt(tem_sigma**2 + tab_sigma**2)

		## COMPUTE EXP(-TAU), TAU & 1-SIGMA STANDARD SEVIATION OF TAU ##
		etaud      = td/tc1665
		taud       = -np.log(etaud)
		tau_sigma  = get_tb_sigma(xd, taud, avmin1, avmax1, avmin2, avmax2)
		etau_sigma = np.abs(etaud)*np.sqrt( (tab_sigma/td)**2 + (tc_er/tc1665)**2 )
		etau_sigma = get_tb_sigma(xd, etaud, avmin1, avmax1, avmin2, avmax2)
		
		# VELOCITY-RANGE & INDEXES #
		xmax_id  = get_vel_index(xd, xmin)   # xmax_index
		xmin_id  = get_vel_index(xd, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id           # Total number of bins
		vrange   = [xmin_id, xmax_id]
		dv       = (xmax-xmin)/num_chnl

		## (FOR FUN) FIT ABSORPTION LINE FOR TAU, V0 and WIDTH ##
		guesspar,base_range   = peak_info(src,pfl)
		lguess                = [tc1665] + guesspar
		x,etaufit,etau,\
		abp,abper,npar,\
		parbase,pname,parinfo = ab_fit(src,xd,td,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2)

		## Tau and Width to cal. N(OH) ##
		tau_fit = []
		v0_fit  = []
		wid_fit = []
		for i in range(1,len(abp),3):
			tau_fit.append(abp[i])
			v0_fit.append(abp[i+1])
			wid_fit.append(abp[i+2])

		## CALCULATE Tex, CHOOSE etaufit OR etau ##
		t_on     = tde + td # On-source Spectrum
		t_on     = t_on[xmin_id:xmax_id]
		xde      = xd[xmin_id:xmax_id]
		e_tau    = etaud[xmin_id:xmax_id]		
		taud     = taud[xmin_id:xmax_id]

		ton      = []
		xe       = []
		etaue    = []
		etauf    = []
		tau      = []
		for i in range(len(t_on)):
			if(taud[i] > 2.*tau_sigma):
				tau.append(taud[i])
				ton.append(t_on[i])
				xe.append(xde[i])
				etaue.append(e_tau[i])
				etauf.append(etaufit[i])

		ton     = np.asarray(ton,   dtype=np.float64)
		xe      = np.asarray(xe,    dtype=np.float64)
		etaue   = np.asarray(etaue, dtype=np.float64)
		tau     = np.asarray(tau,   dtype=np.float64)
		tex     = (ton-trx-(tbg1665 + tc1665)*etaue)/(1.-etaue)

		# Find the vel-range of the peak #
		peak        = get_peak_vel_range(abp)
		npeak       = len(peak)/2	
		tex_peak    = [0.]*npeak
		texsig_peak = [0.]*npeak
		noh_peak    = [0.]*npeak
		stausig     = [0.]*npeak
		noh_sig     = [0.]*npeak

		ltex        = {}
		ltau        = {}
		lstau_sig   = {}
		for j in range(npeak):
			ltex[j]      = []
			ltau[j]      = []
			lstau_sig[j] = []

		# Cal. Tex for each peak #
		for i in range(0, len(xe)):
			for k in range(0,len(peak),2):
				vmin = peak[0+k]
				vmax = peak[1+k]
				if ( (xe[i]>=vmin) and (xe[i]<=vmax) ) :
					ltex[k/2].append(tex[i])
					ltau[k/2].append(tau[i])
					lstau_sig[k/2].append(tau_sigma**2)

		for k in range(npeak):
			if( (len(ltex[k])>0 ) and (np.sum(ltex[k]) > 0.) ):
				tex_peak[k]    = np.mean(ltex[k])
				texsig_peak[k] = np.std(ltex[k])
				noh_peak[k]    = cst*tex_peak[k]*np.sum(ltau[k])*dv
				noh_sig[k]     = noh_peak[k]*np.sqrt( (texsig_peak[k]/tex_peak[k])**2 + np.sum(lstau_sig[k])*dv**2/(np.sum(ltau[k]))**2  )
			else:
				tex_peak[k]    = 0.
				noh_peak[k]    = 0.
				texsig_peak[k] = 0.
				noh_sig[k]     = 0.

		for k in range(npeak):
			if(n<79):
				print '{0}\t{1}\t\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'\
					.format(n,src,round(tau_fit[k],10),round(v0_fit[k],10),round(wid_fit[k],10),round(tex_peak[k],10),round(texsig_peak[k],10),round(noh_peak[k],10),round(noh_sig[k],10) )
			else:
				print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'\
					.format(n,src,round(tau_fit[k],10),round(v0_fit[k],10),round(wid_fit[k],10),round(tex_peak[k],10),round(texsig_peak[k],10),round(noh_peak[k],10),round(noh_sig[k],10) )

## Compute Tex for 1665 line ##
 #
 # params dict data     Data
 # params dict inf408   Info about the Tb_background at 408MHz
 # params int  bd       OH665 or OH667
 #
 # return Tex and N(OH)
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def cal_tex(data,inf408,bd=1):
	bg408    = read_tbg408_healpy()
 	src_list = list(data.la.srcname)

 	for src in src_list:
	 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = read_bins_to_cal_bg(n)
	 	xmin, xmax                                                        = vel_range(n)

	 	# if (src != '3C75'):
			# continue
		if (xmin == 0. and xmax == 0.):
			continue

	 	## Absorption and Emission data ##
	 	xd  = vlsr1
		td  = ab_avg1
		tde = em_avg1
		cst = 3.99757843817
		frq = 1665.402
		pfl = '../data/gauss_1665_peaks.txt'
		if(bd == 2):
			xd  = vlsr2
			td  = ab_avg2
			tde = em_avg2
			cst = 2.21841824609
			frq = 1667.359
			pfl = '../data/gauss_1667_peaks.txt'

		# xd,td = bin_up(xd,td,nbin=1)

		## Background ##
		tbg1665 = 2.8+get_tb_408(ell,bee,inf408.tb_408)*(408./frq)**2.8 # Tbg from 408MHz
		tbg1665 = 2.8+bg408[src]*(408./frq)**2.8 # Tbg from 408MHz
		if(src=='3C123'):
			tbg1665 = 26.

		## BACKGROUND and THEIR UNCERTAINTIES ##
		tc1665,tc_er        = baseline_from_linear_fit(xd, td, avmin1, avmax1, avmin2, avmax2,fit=False)
		bg_off1665,bgoff_er = baseline_from_linear_fit(xd, tde, evmin1, evmax1, evmin2, evmax2,fit=False)
		trx                 = bg_off1665 - tbg1665

		## 1-SIGMA STANDARD SDEVIATION OF Tabspt and Temmission ##
		tab_sigma  = tc_er
		tem_sigma  = bgoff_er
		trx_sigma  = bgoff_er
		ton_sig    = np.sqrt(tem_sigma**2 + tab_sigma**2)

		## COMPUTE EXP(-TAU), TAU & 1-SIGMA STANDARD SEVIATION OF TAU ##
		etaud      = td/tc1665
		taud       = -np.log(etaud)
		tau_sigma  = get_tb_sigma(xd, taud, avmin1, avmax1, avmin2, avmax2)
		etau_sigma = np.abs(etaud)*np.sqrt( (tab_sigma/td)**2 + (tc_er/tc1665)**2 )
		etau_sigma = get_tb_sigma(xd, etaud, avmin1, avmax1, avmin2, avmax2)

		print '**************************'
	 	print '1) Source: '
	 	print '    ', src
	 	print '2) Tcont:'
	 	print '    ',tc1665
	 	print '3) Background of OFF-SOURCE spectrum:'
	 	print '    ',bg_off1665
	 	print '4) Radio contimuum obtained from 408MHz:'
	 	print '    ',tbg1665
	 	print '5) Receiver Temperature, Trx = bg_off1665 - tbg1665:'
	 	print '    ',trx
		
		# VELOCITY-RANGE & INDEXES #
		xmax_id  = get_vel_index(xd, xmin)   # xmax_index
		xmin_id  = get_vel_index(xd, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id           # Total number of bins
		vrange   = [xmin_id, xmax_id]
		dv       = (xmax-xmin)/num_chnl
		print dv
		print dv*frq/300000.
		print dv*2048*frq/300000.

		## (FOR FUN) FIT ABSORPTION LINE FOR TAU, V0 and WIDTH ##
		guesspar,base_range = peak_info(src,fname='../data/gauss_1665_peaks.txt')
		lguess              = [tc1665] + guesspar
		x,etaufit,etau,\
		abp,abper,npar,\
		parbase,pname,parinfo = ab_fit(src,xd,td,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2)

		## PLOT ##
		colors = ['m','g','b','y','c','r','purple','b']
		plt.plot(x,etau, 'b.-', label='data', ms=10)
		plt.plot(x,etaufit,'r-', label='Gaussian fit', lw=2)
		for i in range(2,len(abp),3):
			# plt.axvline(abp[i]-abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color=colors[(i-3)/4])
			# plt.axvline(abp[i]+abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color=colors[(i-3)/4])
			plt.axvline(abp[i]-abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			plt.axvline(abp[i]+abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
		plt.title(src,fontsize=35)
		plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
		plt.ylabel(r'$e^{-\tau}$',fontsize=35)
		plt.legend(loc='upper right')
		plt.grid()
		plt.show()

		## Tau and Width to cal. N(OH) ##
		tau_fit = []
		v0_fit  = []
		wid_fit = []
		for i in range(1,len(abp),3):
			tau_fit.append(abp[i])
			v0_fit.append(abp[i+1])
			wid_fit.append(abp[i+2])

		## CALCULATE Tex, CHOOSE etaufit OR etau ##
		t_on     = tde + td # On-source Spectrum
		t_on     = t_on[xmin_id:xmax_id]
		xde      = xd[xmin_id:xmax_id]
		e_tau    = etaud[xmin_id:xmax_id]		
		taud     = taud[xmin_id:xmax_id]

		ton      = []
		xe       = []
		etaue    = []
		etauf    = []
		tau      = []
		for i in range(len(t_on)):
			if(taud[i] > 3.*tau_sigma):
				tau.append(taud[i])
				ton.append(t_on[i])
				xe.append(xde[i])
				etaue.append(e_tau[i])
				etauf.append(etaufit[i])

		ton     = np.asarray(ton,   dtype=np.float64)
		xe      = np.asarray(xe,    dtype=np.float64)
		etaue   = np.asarray(etaue, dtype=np.float64)
		tau     = np.asarray(tau,   dtype=np.float64)
		tex     = (ton-trx-(tbg1665 + tc1665)*etaue)/(1.-etaue)

		# Find the vel-range of the peak #
		peak        = get_peak_vel_range(abp)
		npeak       = len(peak)/2	
		tex_peak    = [0.]*npeak
		texsig_peak = [0.]*npeak
		noh_peak    = [0.]*npeak
		stausig     = [0.]*npeak
		noh_sig     = [0.]*npeak

		ltex        = {}
		ltau        = {}
		lstau_sig   = {}
		for j in range(npeak):
			ltex[j]      = []
			ltau[j]      = []
			lstau_sig[j] = []

		# Cal. Tex for each peak #
		for i in range(0, len(xe)):
			for k in range(0,len(peak),2):
				vmin = peak[0+k]
				vmax = peak[1+k]
				if ( (xe[i]>=vmin) and (xe[i]<=vmax) ) :
					ltex[k/2].append(tex[i])
					ltau[k/2].append(tau[i])
					lstau_sig[k/2].append(tau_sigma**2)

		for k in range(npeak):
			if( (len(ltex[k])>0 ) and (np.sum(ltex[k]) > 0.) ):
				tex_peak[k]    = np.mean(ltex[k])
				texsig_peak[k] = np.std(ltex[k])
				noh_peak[k]    = cst*tex_peak[k]*np.sum(ltau[k])*dv
				noh_sig[k]     = noh_peak[k]*np.sqrt( (texsig_peak[k]/tex_peak[k])**2 + np.sum(lstau_sig[k])*dv**2/(np.sum(ltau[k]))**2  )
			else:
				tex_peak[k]    = 0.
				noh_peak[k]    = 0.
				texsig_peak[k] = 0.
				noh_sig[k]     = 0.

		print ''
		print '*** Tex & N(OH) of each peak: Tex, Tex_er, NOH, NOH_er ***'
		snoh = 0.
		for k in range(npeak):
			print 'Peak ' + str(k) + ' :'
			print '    ',round(tex_peak[k],10),round(texsig_peak[k],10),round(noh_peak[k],10),round(noh_sig[k],10)

		## PLOT On-Source LINE ##
		plt.plot(xe,ton)
		for i in range(0,len(peak),2):
			plt.axvline(peak[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			plt.axvline(peak[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
		plt.title(src + ' - On-Source line')
		plt.grid()
		plt.show()

		## PLOT ##
		plt.plot(xe,tex, 'r', label='$T_{ex}$', lw=2)
		plt.plot(xe,-7000.*np.log(etaue), 'b', label=r'$7000*e^{-\tau}$', lw=2)
		for i in range(0,len(peak),2):
		# 	plt.axvline(peak[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
		# 	plt.axvline(peak[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			plt.axvline(peak[i], ymin=-10., ymax=1000.,linewidth=3,color='k')
			plt.axvline(peak[i+1], ymin=-10., ymax=1000.,linewidth=3,color='k')
		# for i in range(2,len(abp),3):
		# 	plt.axvline(abp[i]-abp[i+1]/4./np.sqrt(2),ymin=-10., ymax=1000.,linewidth=1,color='k') # Sigma range of Gaussian
		# 	plt.axvline(abp[i]+abp[i+1]/4./np.sqrt(2),ymin=-10., ymax=1000.,linewidth=1,color='k')
		plt.ylim(-50.,400.)
		plt.title(src, fontsize=35 )
		plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
		plt.ylabel('$T_{ex} (K)$',fontsize=35)
		plt.legend(loc='upper right')
		plt.grid()
		plt.show()

	    ## CALCULATE N(OH) ##
		plt.plot(xe,-np.log(etaue))
		for i in range(0,len(peak),2):
			plt.axvline(peak[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			plt.axvline(peak[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])

		plt.title(src + ' - Opacity, tau')
		plt.grid()
		plt.show()

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408, Continuum at 408MHz
# cal_tex(data, inf408, bd=1)
cal_tex_print(data, inf408, bd=1)

sys.exit()