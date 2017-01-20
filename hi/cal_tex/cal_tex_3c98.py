import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot
from gauss_fit      import gfit

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

## Calculate multiple (N) Gaussians + offset. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def gcurv(xdata, zro1, hgt1, cen1, wid1):
	#DETERMINE NR OF GAUSSIANS...
	ngaussians = 1 # cho vui

	tfit = 0.*xdata + zro1
	if (wid1 > 0.):
		tfit = tfit + hgt1*np.exp(- ( (xdata-cen1)/(0.6005612*wid1))**2)

	return tfit

## Convert angle to a specified range by finding the angle modulo the extent of the range. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def tb_exp(xdata, zrocnm, hgtcnm, cencnm, widcnm, tspincnm, ordercnm, continuum, hgtwnm, cenwnm, widwnm, fwnm):
	#ZEROTH STEP IS TO REARRANGE CLOUDS IN ORDER OF 'ORDER'.
	zro1   = zrocnm
	hgt1   = list(hgtcnm)
	cen1   = list(cencnm)
	wid1   = list(widcnm)
	tspin1 = list(tspincnm)

	#FIRST STEP IS TO CALCULATE THE OPACITY OF EACH COLD CLOUD...
	nrcnm  = len(hgt1)
	taucnm = np.zeros((len(xdata), nrcnm))

	for nrc in range(nrcnm):
		tau1nrc = gcurv(xdata, zro1, hgt1[nrc], cen1[nrc], wid1[nrc])
		taucnm[:, nrc] = tau1nrc

	if (len(ordercnm) != 1):
		tausum = taucnm.sum(1)
	else:
		tausum = taucnm

	exp_tausum = np.exp(-tausum)

	##********** CALCULATE THE WNM CONTRIBUTION ********************
	##  EXPRESS THE WNM CONTRIBUTION AS A SUM OF GAUSSIANS:
	##	FWNM, ZROWNM, HGTWNM, CENWNM, WIDWNM
	tb_cont    = continuum* exp_tausum

	tb_wnm_tot = np.zeros(len(xdata))
	nrwnm      = len(hgtwnm)
	for nrw in range(nrwnm):
		tb_wnm_nrw = gcurv(xdata, 0., hgtwnm[nrw], cenwnm[nrw], widwnm[nrw])
		tb_wnm_tot = tb_wnm_tot + tb_wnm_nrw*(fwnm[nrw] + (1.-fwnm[nrw])*exp_tausum)

	#*************** NEXT CALCULATE THE CNM CONTRIBUTION ****************

	tb_cnm_tot = np.zeros(len(xdata))

	# BRIGHT TEMP OF EACH CNM CLUMP:
	tbclump = np.zeros((len(xdata), nrcnm))
	for nrc in range(nrcnm):
		tbclump[:,nrc] = tspin1[nrc] * (1. - np.exp(-taucnm[:,nrc]))

	for nrc in range(nrcnm):
		temp        = np.reshape(taucnm[:, 0:nrc+1], (len(xdata), nrc+1))
		tausum_nrc  = temp.sum(1)
		exp_tau_nrc = np.exp(taucnm[:, nrc] - tausum_nrc)
		tb_cnm_tot  = tb_cnm_tot + tspin1[nrc] * (1. - np.exp(-taucnm[:,nrc]) ) * exp_tau_nrc

	tb_tot = tb_cont+ tb_cnm_tot + tb_wnm_tot

	return tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tausum

## Get genaral info of the HI source ##
 #
 # params dict data Data
 # params str src Source-name
 #
 # return general Infos of source
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_src_hi_info(data, src, src_list):
	n      = src_list.index(src)

	ra50   = data.la.ra1950
	dec50  = data.la.dec1950
	ell    = data.la.ell
	bee    = data.la.bee

	freq   = data.la.cfr_bd0
	vlsr   = data.la.vlsr_bd0

	em_avg = 0.5*correct_ctrl_chnl(data.la.i_em_avg_bd0)
	em_med = 0.5*correct_ctrl_chnl(data.la.i_em_med_bd0)
	ab_avg = 0.5*correct_ctrl_chnl(data.la.i_abs_avg_bd0)
	ab_med = 0.5*correct_ctrl_chnl(data.la.i_abs_med_bd0)

	return n,ell[n],bee[n],freq[n],vlsr[n], em_avg[n],ab_avg[n]

## Read vel-range to calculate background ##
 #
 # params int n Order of the Source
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_bins_to_cal_bg(n,fname='../data/bins_to_cal_bg_for_hi.txt'):
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

## Velocity range ##
 #
 # params int n order of Source
 # params str fname Filename
 #
 # return float xmin, xmax Vel-range to fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def vel_range(n, fname = '../data/vel_range_ms101src.txt'):
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']	
	
	xmin = vmin[n]
	xmax = vmax[n]

	return xmin, xmax

## Read peaks Info, initial params for fitting ##
 #
 # params string fname File-name
 # params string sc    Source

 # return list guessp Guess-Parameters 
 #        list base_range Range of Baselines
 #        float baseline Baselines after checking chi2
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def peak_guessp(sc,fname='../data/gauss_peaks.txt'):
	ret = {}

	cols = ['idx','src','tau','v0','wid','bmin','bmax', 'v1', 'v2', 'tc', 'bgoff', 'tbg', 'bsl']
	fmt  = ['i',  's',   'f', 'f', 'f',   'f',   'f',    'f',   'f',    'f',   'f',  'f',   'f']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read()
	bmin = info['bmin']
	bmax = info['bmax']
	bsl  = info['bsl']
	src  = info['src']
	tau  = info['tau']
	v0   = info['v0']
	wid  = info['wid']
	vid1 = info['v1']
	vid2 = info['v2']
	tc   = info['tc']
	bge  = info['bgoff']
	tbg  = info['tbg']
	for i in range(0,len(src)):
		if info['src'][i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['base_range'] = [bmin[i],bmax[i]]
			ret[src[i]]['baseline']   = bsl[i]
			ret[src[i]]['tau']        = [tau[i]]
			ret[src[i]]['v0']         = [v0[i]]
			ret[src[i]]['wid']        = [wid[i]]
			ret[src[i]]['vid1']       = vid1[i]
			ret[src[i]]['vid2']       = vid2[i]
			ret[src[i]]['tc']         = tc[i]
			ret[src[i]]['bge']        = bge[i]
			ret[src[i]]['tbg']        = tbg[i]
		else:
			ret[src[i]]['base_range'] = [bmin[i],bmax[i]]
			ret[src[i]]['baseline']   = bsl[i]
			ret[src[i]]['tau']        += [tau[i]]
			ret[src[i]]['v0']         += [v0[i]]
			ret[src[i]]['wid']        += [wid[i]]
			ret[src[i]]['vid1']       = vid1[i]
			ret[src[i]]['vid2']       = vid2[i]
			ret[src[i]]['tc']         = tc[i]
			ret[src[i]]['bge']        = bge[i]
			ret[src[i]]['tbg']        = tbg[i]

	return ret[sc]['tau'], ret[sc]['v0'], ret[sc]['wid'], ret[sc]['base_range'], \
			ret[sc]['vid1'], ret[sc]['vid2'], ret[sc]['tc'], ret[sc]['bge'], ret[sc]['tbg'], ret[sc]['baseline']

## Read vel-range to calculate background ##
 #
 # params 1-D array x x-data
 # params 1-D array y T-data
 # params float vmin1 Range of velocity
 # params float vmax1 Range of velocity
 # params float vmin2 Range of velocity
 # params float vmax2 Range of velocity
 #
 # return slope and intercept of the Linear fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def baseline_from_linear_fit(x,y,vmin1,vmax1,vmin2,vmax2,fit=True):
	if(fit):
		tb = []
		v  = []
		count = 0
		for i in range(0,len(x)):
			if (((x[i] > vmin1) and (x[i] < vmax1)) or ((x[i] > vmin2) and (x[i] < vmax2))):
				tb.append(y[i])
				v.append(x[i])
				count = count + 1

	 	slope,bsline = np.polyfit(v, tb, 1)
	else:
	 	bsline = 0.
		count  = 0.
		slope  = 0.
		for i in range(0,len(x)):
			if (((x[i] > vmin1) and (x[i] < vmax1)) or ((x[i] > vmin2) and (x[i] < vmax2))):
				bsline = bsline + y[i]
				count  = count + 1

		bsline = bsline/count

 	sigma = 0.
 	for i in range(0,len(x)):
		if (((x[i] > vmin1) and (x[i] < vmax1)) or ((x[i] > vmin2) and (x[i] < vmax2))):
			sigma = sigma + (y[i]-(slope*x[i] + bsline))**2

	sigma = np.sqrt(sigma/(count-1))
 	return bsline, sigma

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
	fit      = gfit()
 	src_list = list(data.la.srcname)
 	bg408    = read_tbg408_healpy('../../oh/result/bg408_to_compare.txt')

	for src in src_list:
		if(src != '3C98'):
			continue

		n = src_list.index(src)
		n,ell,bee,freq,vlsr, em_avg,ab_avg = get_src_hi_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,\
	 	evmin1,evmax1,evmin2,evmax2        = read_bins_to_cal_bg(n, '../data/bins_to_cal_bg_for_hi.txt')
	 	xmin, xmax                         = vel_range(n)
		if (xmin == 0. and xmax == 0.):
			continue

		print n, '   ', src, ell, bee, freq

		v  = vlsr
		ta = ab_avg
		te = em_avg

		## COntinuum Background
		# tbg = 2.7+bg408[src]*(408./1420.405)**2.7 # Tbg from 408MHz, chua co so lieu cho 3c245
		tbg = 2.7+get_tb_408(ell,bee,inf408.tb_408)*(408./1420.405)**2.7 # Tbg from 408MHz

		## BACKGROUNDS and THEIR UNCERTAINTIES ##
		tc,tc_er        = baseline_from_linear_fit(v, ta, avmin1, avmax1, avmin2, avmax2, fit=False)
		bg_off,bgoff_er = baseline_from_linear_fit(v, te, evmin1, evmax1, evmin2, evmax2, fit=False)
		tau, v0, wid, base_range, v1,v2, tc1665, bg_off1665, tbg1665, base = peak_guessp(src,'../data/gauss_peaks.txt')

		print n, '   ', src, ell, bee
		print 'Tc: ', tc, tc_er
		print 'Tbg_off: ', bg_off, bgoff_er
		print 'Tbg: ', tbg

		# Get Index of the Velocity-range #
		xmax_id  = get_vel_index(v, xmin)   # xmax_index
		xmin_id  = get_vel_index(v, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id          # Total number of bins
		vrange   = [xmin_id, xmax_id]
		print 'vrange', vrange

		print n, '   ', src, ell, bee, tbg
		continuum_em = tbg

		look   = -1
		zro0yn = 1
		hgt0yn = [1]
		cen0yn = [1]
		wid0yn = [1]
		npeak  = 4

		continuum_defln = tc #28.29
		continuumyn = 1

		zrocnm   = 0.
		taucnm   = [0.51, 0.297, 0.041, 0.092]
		cencnm   = [9.5,9.5,23.0,-1.0]
		widcnm   = [5.21, 1.35, 6.82, 4.87]
		tspincnm = [1.e-6]*npeak
		ordercnm = list(range(npeak))

		zrocnmyn   = 0
		taucnmyn   = [1]*npeak
		cencnmyn   = [1]*npeak
		widcnmyn   = [1]*npeak
		tspincnmyn = [0]*npeak

		hgtwnm = [0.003]
		cenwnm = [-10.4]
		widwnm = [21.46]
		fwnm   = [0.5]

		continuum_deflnyn = 1
		hgtwnmyn = 0
		cenwnmyn = 0
		widwnmyn = 0
		fwnmyn   = 0

		# BD1 ===== 1665  bd1
		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(v,zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
												continuum_defln, hgtwnm, cenwnm, widwnm, fwnm)

		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_defln1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, v, ta, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
        	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
        	continuum_defln, hgtwnm, cenwnm, widwnm, fwnm, \
        	continuum_deflnyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(v, \
			zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
			continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		nrgauss_cnm = len(taucnm1)
		nrg_cnm     = nrgauss_cnm

		## Absoprtion line
		# plt.plot(v,ta, 'b-', linewidth=1, label='data, Absorption line')
		# plt.plot(v,tb_cont, 'r-', linewidth=1, label='data, Absorption line')
		# plt.title(src, fontsize=30)
		# plt.ylabel('$T_{b} [K]$', fontsize=35)
		# plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
		# # plt.xlim(0.0, 2.0)
		# # plt.xlim(-1.0, 6.0)
		# plt.grid(True)
		# plt.tick_params(axis='x', labelsize=18)
		# plt.tick_params(axis='y', labelsize=15)

		# # plt.text(0.0, 3.2, r'$f = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al.', color='blue', fontsize=20)
		# plt.legend(loc='upper left', fontsize=18)
		# plt.show()

		# ======== Do emission spectrum =================                                         
		tspincnm   = [97.27, 27.59, 134.11, 16.82]
		tspincnmyn = [1]*npeak
		tbaseline  = bg_off

		# Tst 7 #
		zrocnm = zrocnm1
		taucnm = taucnm1
		cencnm = cencnm1
		widcnm = widcnm1
		nrcnm  = len(cencnm1)

		zrocnmyn       = 0
		continuum_emyn = 0
		taucnmyn       = [0]*nrcnm
		cencnmyn       = [0]*nrcnm
		widcnmyn       = [0]*nrcnm
		ordercnm       = list(range(nrcnm))

		## WNM
		hgtwnm = [0.021]
		cenwnm = [6.3]
		widwnm = [14.86]
		fwnm   = [0.53]

		nrwnm    = len(hgtwnm)
		hgtwnmyn = [1]*nrwnm
		cenwnmyn = [1]*nrwnm
		widwnmyn = [1]*nrwnm
		fwnmyn   = [0]*nrwnm

		tee = te - tbaseline + tbg
		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(v, \
				zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
				continuum_em, hgtwnm, cenwnm, widwnm, fwnm)

		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_em1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, v, tee, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
	    	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
	    	continuum_em, hgtwnm, cenwnm, widwnm, fwnm, \
	    	continuum_emyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(v, \
				zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
				continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		nrgauss_wnm = len(hgtwnm1)
		nrg_wnm     = nrgauss_wnm
		# End  - Tst 7 #

		print 'bd2...problem: ', problem

		print '####### HI ######'
		print '1) Tex & Error:'
		print '    ', tspincnm1
		print '    ', sigtspincnm1
		print '2) Tau & Error:'
		print '    ', taucnm1
		print '    ', sigtaucnm1
		print '3) V0 & Error:'
		print '    ', cencnm1
		print '    ', sigcencnm1
		print '4) Width & Error:'
		print '    ', widcnm1
		print '    ', sigwidcnm1
		print '5) hgtwnm & Error:'
		print '    ', hgtwnm1
		print '    ', sighgtwnm1
		print '6) cenwnm & Error:'
		print '    ', cenwnm1
		print '    ', sigcenwnm1
		print '7) widwnm & Error:'
		print '    ', widwnm1
		print '    ', sigwidwnm1
		print '8) fwnm & Error:'
		print '    ', fwnm1
		print '    ', sigfwnm1
		print '9) nloop:'
		print '    ', nloop
#============== MAIN ==============#
data   = readsav('../../oh/data/makelines.sav') #data.la
inf408 = readsav('../../oh/data/tb_408.sav') # l_cntr, b_cntr, tb_408
cal_tex(data, inf408)

sys.exit()