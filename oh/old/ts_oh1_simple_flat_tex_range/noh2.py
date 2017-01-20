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
	cols     = ['src','l','b', 'il', 'ib', 'tbg', 'l-idx','b-idx','tbg1','tbg_hpy']
	fmt      = ['s','f','f','f','f','f','f','f','f','f']
	src      = restore(fname, 2, cols, fmt)
	info     = src.read()

	src = info['src']
	gl  = info['l']
	gb  = info['b']
	il  = info['il']
	ib  = info['ib']
	tbg = info['tbg']
	lid = info['l-idx']
	bid = info['b-idx']
	bg1 = info['tbg1']
	bgh = info['tbg_hpy']

	ret = {}
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
	for i in range(len(xd)):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			tb.append(td[i])
			v.append(xd[i])

	## Linear fit ##
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

	## CALCULATE 1-SIGMA UNCERTAINTY OF Tb DATA, STANDARD DEVIATION ##
	sigma = 0.
	n     = len(tbe)
	for i in range(n):
		sigma = sigma + ((tbe[i] - u)**2)/(n-1)

 	sigma = np.sqrt(sigma)

 	ery = np.zeros(td.shape)
 	for i in range(len(xd)):
		ery[i] = sigma

 	return sigma,ery

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

	ret = {}
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
def get_src_info(data,src,src_list):
	n     = src_list.index(src)

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

## Compute the Tex spectrum from etau ##
 #
 # params list etau      Exp(-tau)
 # params list velo      List of Velocity
 # params list tb_off    List of Off-source Temperature
 # params float tbg      CMB + Radio noise
 # params float trx      Receiver Temperature
 # params list popt      Fit-result parameters 
 # params 1-D Array tau  Tau
 # params float dv       Velocity-bin
 # params float cst      Factor to calculate N(OH) from Tex
 # params 1-D Array tex_sig  Uncertainties of Tex
 # params 1-D Array tau_sig  Uncertainties of Tau
 #
 # return List tex       Tex spectrum
 # return float s        Tex Mean
 # return int count      Number of channels to compute the Tex
 # return list tex_peak  Tex of each peak
 # return list tex_count Number of channels to compute the Tex for each peak
 # return list peak      Vel-range of peaks at FWHM
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tex_from_etau(etau,velo,tb_off,tbg,trx,popt,tau,dv,cst,tex_sig,tausig):
	# Compute the Excitation Temperature #
	tex   = []
	tboff = [] 
	for i in range(0, len(etau)):
		tboff.append(tb_off[i]-trx)
		if (etau[i] == 1.):
			tex.append(0.)
		else:
			ts = (tb_off[i]-trx-tbg*etau[i])/(1.-etau[i])
			tex.append(ts)
	
	# Find the vel-range of the peak #
	peak      = get_peak_vel_range(popt)
	npeak     = len(peak)/2
	tex_peak  = [0.]*npeak
	stau_peak = [0.]*npeak
	tex_count = [0]*npeak
	ts_sig    = [0.]*npeak
	stausig   = [0.]*npeak
	noh_sig   = [0.]*npeak

	# Cal. Tex for each peak #
	s         = 0.
	count     = 0
	stau_sig  = 0.
	for i in range(0, len(velo)):
		for k in range(0,len(peak),2):
			vmin = peak[0+k]
			vmax = peak[1+k]
			if ( (velo[i]>=vmin) and (velo[i]<=vmax) ) :
				tex_peak[k/2]  = tex_peak[k/2] + tex[i]
				stau_peak[k/2] = stau_peak[k/2] + tau[i]*dv
				tex_count[k/2] = tex_count[k/2] + 1
				ts_sig[k/2]    = ts_sig[k/2] + tex_sig[i]
				stausig[k/2]   = stausig[k/2] + tausig[i]**2
				s              = s + tex[i]
				count          = count + 1

	ss = 0.
	for k in range(0,len(peak)/2):
		if(tex_peak[k] > 0.):
			tex_peak[k] = tex_peak[k]/tex_count[k]
			ts_sig[k]   = ts_sig[k]/tex_count[k]
			ds2   = stausig[k] + (stau_peak[k]/dv)**2
			dtex2 = (ts_sig[k]/tex_peak[k])**2
		else:
			tex_peak[k] = 0.
			ts_sig[k]   = 0.
			ds2   = stausig[k] + (stau_peak[k]/dv)**2
			dtex2 = 0.
		
		
		stau_peak[k] = cst*tex_peak[k]*stau_peak[k]
		noh_sig[k]   = stau_peak[k]*np.sqrt(ds2+dtex2)
		ss           = ss + tex_peak[k]
	
	s = s/count
	return tex,s,count,tex_peak,tex_count,peak,tboff,stau_peak,ts_sig,noh_sig

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

	##  data and 1-sigma uncertainty ##
	x    = xd[xmin_id:xmax_id]
	y    = td[xmin_id:xmax_id]
	sg,\
	erry = get_tb_sigma(xd, td, evmin1, evmax1, evmin2, evmax2)
	erry = erry[xmin_id:xmax_id]

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

 	## Carl's Tex ##
	ts_carl  = read_tex_carl()

 	for src in src_list:
	 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = read_bins_to_cal_bg(n)
	 	xmin, xmax                                                        = vel_range(n)

		if (xmin == 0. and xmax == 0.):
			continue

	 	## Absorption and Emission data ##
	 	xd  = vlsr1
		td  = ab_avg1
		tde = em_avg1
		cst = 2.39854792704
		frq = 1665.402
		pfl = '../data/gauss_1665_peaks.txt'
		if(bd == 2):
			xd  = vlsr2
			td  = ab_avg2
			tde = em_avg2
			cst = 2.21841851219
			frq = 1667.359
			pfl = '../data/gauss_1667_peaks.txt'

		## Background ##
		tbg1665    = 2.8+get_tb_408(ell,bee,inf408.tb_408)*(408./frq)**2.8 # Tbg from 408MHz
		tbg1665    = 2.8+bg408[src]*(408./frq)**2.8 # Tbg from 408MHz
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

		## COMPUTE EXP(-TAU), TAU & 1-SIGMA STANDARD SEVIATION OF TAU ##
		etaud      = td/tc1665
		taud       = -np.log(etaud)
		tau_sigma  = get_tb_sigma(xd, taud, avmin1, avmax1, avmin2, avmax2)
		etau_sigma = np.abs(etaud)*np.sqrt( (tab_sigma/td)**2 + (tc_er/tc1665)**2 )
		
		# VELOCITY-RANGE & INDEXES #
		xmax_id  = get_vel_index(xd, xmin)   # xmax_index
		xmin_id  = get_vel_index(xd, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id           # Total number of bins
		vrange   = [xmin_id, xmax_id]
		dv       = (xmax-xmin)/num_chnl

		## FIT ABSORPTION LINE FOR TAU, V0 and WIDTH ##
		guesspar,base_range = peak_info(src,fname=pfl)
		lguess              = [tc1665] + guesspar
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
		tde        = tde[xmin_id:xmax_id]
		xde        = xd[xmin_id:xmax_id]
		e_tau      = etaud[xmin_id:xmax_id]		
		taud       = taud[xmin_id:xmax_id]
		etau_sigma = etau_sigma[xmin_id:xmax_id]
		tau_sig    = np.abs(etau_sigma/e_tau)

		te       = []
		xe       = []
		etaue    = []
		tau      = []
		etau_sig = []
		tausig   = []
		for i in range(len(tde)):
			if(taud[i] > 2.*tau_sigma[0]):
				tau.append(taud[i])
				te.append(tde[i])
				xe.append(xde[i])
				etaue.append(e_tau[i])
				tau.append(taud[i])
				etau_sig.append(etau_sigma[i])
				tausig.append(tau_sig[i])

		te       = np.asarray(te, dtype=np.float64)
		xe       = np.asarray(xe, dtype=np.float64)
		etaue    = np.asarray(etaue, dtype=np.float64)
		tau      = np.asarray(tau, dtype=np.float64)
		etau_sig = np.asarray(etau_sig, dtype=np.float64)
		tausig   = np.asarray(tausig, dtype=np.float64)

		a = np.abs((te-bg_off1665)/(1.-etaue))
		b = np.sqrt(2.*tem_sigma**2)/(te-bg_off1665)
		c = etau_sig/(1.-etaue)
		d = np.sqrt(b**2 + c**2)
		u = a*d

		# plt.plot(xe,u)
		# # plt.plot(xde,u/tex)
		# # plt.plot(xde,u)
		# plt.show()
		# sys.exit()

		tex,s,count,tex_peak,tex_count,peak,\
		tboff,noh_peak,ts_sig, noh_sig = cal_tex_from_etau(etaue,xe,te,tbg1665,trx,abp,tau,dv,cst,u,tausig)



		## Carl's Tex ##
		ctex  = ts_carl[src]
		ctau  = ctex['tau']
		ctsp  = ctex['ts_carl']

		snoh  = 0.
		diff  = 0.
		count = 0
		for k in range(0,len(peak)/2):
			if(n<79):
				print '{0}\t{1}\t\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'\
					.format(n,src,round(tau_fit[k],10),round(v0_fit[k],10),round(wid_fit[k],10),round(tex_peak[k],10),round(ts_sig[k],10),round(noh_peak[k],10),round(noh_sig[k],10) )
			else:
				print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'\
					.format(n,src,round(tau_fit[k],10),round(v0_fit[k],10),round(wid_fit[k],10),round(tex_peak[k],10),round(ts_sig[k],10),round(noh_peak[k],10),round(noh_sig[k],10) )

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

	 	# if (src != '3C454.3'):
			# continue
		if (xmin == 0. and xmax == 0.):
			continue

	 	## Absorption and Emission data ##
	 	xd  = vlsr1
		td  = ab_avg1
		tde = em_avg1
		cst = 2.39854792704
		frq = 1665.402
		pfl = '../data/gauss_1665_peaks.txt'
		if(bd == 2):
			xd  = vlsr2
			td  = ab_avg2
			tde = em_avg2
			cst = 2.21841851219
			frq = 1667.359
			pfl = '../data/gauss_1667_peaks.txt'

		# xd,td = bin_up(xd,td,nbin=1)

		## Background ##
		tbg1665    = 2.8+get_tb_408(ell,bee,inf408.tb_408)*(408./frq)**2.8 # Tbg from 408MHz
		tbg1665    = 2.8+bg408[src]*(408./frq)**2.8 # Tbg from 408MHz
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

		## COMPUTE EXP(-TAU), TAU & 1-SIGMA STANDARD SEVIATION OF TAU ##
		etaud      = td/tc1665
		taud       = -np.log(etaud)
		tau_sigma  = get_tb_sigma(xd, taud, avmin1, avmax1, avmin2, avmax2)
		etau_sigma = np.abs(etaud)*np.sqrt( (tab_sigma/td)**2 + (tc_er/tc1665)**2 )

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

		## (FOR FUN) FIT ABSORPTION LINE FOR TAU, V0 and WIDTH ##
		guesspar,base_range = peak_info(src,fname='../data/gauss_1665_peaks.txt')
		lguess              = [tc1665] + guesspar
		x,etaufit,etau,\
		abp,abper,npar,\
		parbase,pname,parinfo = ab_fit(src,xd,td,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2)

		## PLOT ##
		colors = ['m','g','b','y','c','r','purple','b']
		plt.plot(x,etau)
		plt.plot(x,etaufit,'r-')
		for i in range(2,len(abp),3):
			plt.axvline(abp[i]-abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color=colors[(i-3)/4])
			plt.axvline(abp[i]+abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color=colors[(i-3)/4])
		plt.title('Absorption fit - ' + src)
		plt.grid()
		plt.show()

		## Tau and Width to cal. N(OH) ##
		tau_fit = []
		wid_fit = []
		for i in range(1,len(abp),3):
			tau_fit.append(abp[i])
			wid_fit.append(abp[i+2])

		## CALCULATE Tex, CHOOSE etaufit OR etau ##
		tde        = tde[xmin_id:xmax_id]
		xde        = xd[xmin_id:xmax_id]
		e_tau      = etaud[xmin_id:xmax_id]		
		taud       = taud[xmin_id:xmax_id]
		etau_sigma = etau_sigma[xmin_id:xmax_id]
		tau_sig    = np.abs(etau_sigma/e_tau)

		te       = []
		xe       = []
		etaue    = []
		etauf    = []
		tau      = []
		etau_sig = []
		tausig   = []
		for i in range(len(tde)):
			if(taud[i] > 2.*tau_sigma[0]):
				tau.append(taud[i])
				te.append(tde[i])
				xe.append(xde[i])
				etaue.append(e_tau[i])
				etauf.append(etaufit[i])
				tau.append(taud[i])
				etau_sig.append(etau_sigma[i])
				tausig.append(tau_sig[i])

		te       = np.asarray(te, dtype=np.float64)
		xe       = np.asarray(xe, dtype=np.float64)
		etaue    = np.asarray(etaue, dtype=np.float64)
		tau      = np.asarray(tau, dtype=np.float64)
		etau_sig = np.asarray(etau_sig, dtype=np.float64)
		tausig   = np.asarray(tausig, dtype=np.float64)

		a = np.abs((te-bg_off1665)/(1.-etaue))
		b = np.sqrt(2.*tem_sigma**2)/(te-bg_off1665)
		c = etau_sig/(1.-etaue)
		d = np.sqrt(b**2 + c**2)
		u = a*d

		# plt.plot(xe,u)
		# # plt.plot(xde,u/tex)
		# # plt.plot(xde,u)
		# plt.show()
		# sys.exit()

		tex,s,count,tex_peak,tex_count,peak,\
		tboff,noh_peak,ts_sig, noh_sig = cal_tex_from_etau(etaue,xe,te,tbg1665,trx,abp,tau,dv,cst,u,tausig)

		print ''
		print '*** Tex & N(OH) of each peak: Tex, Tex_count,Tau_fit,NOH ***'
		snoh = 0.
		for k in range(0,len(peak)/2):
			print 'Peak ' + str(k) + ' :'
			print '    ',tex_peak[k], tex_count[k],tau_fit[k],noh_peak[k]
			snoh = snoh + noh_peak[k]

		# print '=> Mean Tex:'
		# print '    ',s

		## PLOT EMISSION LINE ##
		plt.plot(xe,tboff)
		for i in range(0,len(peak),2):
			plt.axvline(peak[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			plt.axvline(peak[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
		plt.title(src + ' - Emission line')
		plt.grid()
		plt.show()

		## PLOT ##
		plt.plot(xe,tex)
		plt.plot(xe,-7000.*np.log(etaue))
		for i in range(0,len(peak),2):
			plt.axvline(peak[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			plt.axvline(peak[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
		for i in range(2,len(abp),3):
			plt.axvline(abp[i]-abp[i+1]/4./np.sqrt(2),ymin=-10., ymax=1000.,linewidth=1,color='k') # Sigma range of Gaussian
			plt.axvline(abp[i]+abp[i+1]/4./np.sqrt(2),ymin=-10., ymax=1000.,linewidth=1,color='k')
		plt.ylim(-50.,400.)
		plt.title(src + ' - Tex')
		plt.grid()
		plt.show()

	    ## CALCULATE N(OH) ##
		plt.plot(xe,-np.log(etaue))
		plt.plot(xe,-np.log(etauf))
		for i in range(0,len(peak),2):
			plt.axvline(peak[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			plt.axvline(peak[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])

		plt.title(src + ' - Opacity, tau')
		plt.grid()
		plt.show()

		popt = []
		for i in range(1, len(abp), 3):
			popt += [abp[i], abp[i+1], abp[i+2], tex_peak[(i-1)/3]]

		stau_fit = quad(tau_func, xmin, xmax, args=tuple(popt))
		print '0) Tex[i]*integral(tau)'
		print '    ', stau_fit

		noh_fit = cst*stau_fit[0]  # x10^14
		print '1) N(OH):'
		print '    ', noh_fit, 'x10^14'
		print '2) Sum of N(OH):'
		print '    ',snoh, 'x10^14'

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408
# cal_tex(data, inf408, bd=2)
cal_tex_print(data, inf408, bd=2)

sys.exit()