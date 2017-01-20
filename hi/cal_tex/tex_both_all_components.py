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

	# em_avg1 = 0.5*correct_ctrl_chnl(data.la.i_em_avg_bd1)
	# em_med1 = 0.5*correct_ctrl_chnl(data.la.i_em_med_bd1)
	# ab_avg1 = 0.5*correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	# ab_med1 = 0.5*correct_ctrl_chnl(data.la.i_abs_med_bd1)

	# em_avg2 = 0.5*correct_ctrl_chnl(data.la.i_em_avg_bd2)
	# em_med2 = 0.5*correct_ctrl_chnl(data.la.i_em_med_bd2)
	# ab_avg2 = 0.5*correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	# ab_med2 = 0.5*correct_ctrl_chnl(data.la.i_abs_med_bd2)

	em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	return n,ell[n],bee[n],oh_f1[n],vlsr1[n],oh_f2[n],vlsr2[n],em_avg1[n],ab_avg1[n],em_avg2[n],ab_avg2[n]

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

## Velocity range ##
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

	# Gaussian peaks' info - estimate values of Amp, tau0 and Width #
	peak = peak_info('../data/gauss_1665_peaks.txt')

	for src in src_list:
		if(src != '3C123'):
			continue

		n = src_list.index(src)
		n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = read_bins_to_cal_bg(n)
	 	xmin, xmax                                                        = vel_range(n)

		if (xmin == 0. and xmax == 0.):
			continue

		tbg=50./2
		continuum_em= tbg

		# VLSR #
		x = vlsr1
		velplotmin = xmin
		velplotmax = xmax

		# Get Index of the Velocity-range #
		xmax_id  = get_vel_index(x, xmin)   # xmax_index
		xmin_id  = get_vel_index(x, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id          # Total number of bins
		vrange   = [xmin_id, xmax_id]
		vrange   = [966,1023]

		print n, '   ', src, ell, bee

		look   = -1
		zro0yn = 1
		hgt0yn = [1,1,1]
		cen0yn = [1,1,1]
		wid0yn = [1,1,1]
		npeak  = 3

		continuum_defln = 532.
		continuumyn = 1

		zrocnm   = 0.
		taucnm   = [30./continuum_defln, 30./continuum_defln, 30./continuum_defln]
		cencnm   = [3.2, 4.4, 5.3]
		widcnm   = [0.6,0.6,0.6]
		tspincnm = [1.e-6]*npeak
		ordercnm = list(range(npeak))

		zrocnmyn   = 0
		taucnmyn   = [1]*npeak
		cencnmyn   = [1]*npeak
		widcnmyn   = [1]*npeak
		tspincnmyn = [0]*npeak

		hgtwnm = [1.0e-6]
		cenwnm = [30.]
		widwnm = [4.]
		fwnm   = [1.0]

		continuum_deflnyn = 1
		hgtwnmyn = 0
		cenwnmyn = 0
		widwnmyn = 0
		fwnmyn   = 0

		# BD1 ===== 1665  bd1
		x = vlsr1
		t = ab_avg1

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x,zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
												continuum_defln, hgtwnm, cenwnm, widwnm, fwnm)

		plt.plot(x,t)
		plt.plot(x, tb_cont, 'g-')

		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_defln1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, x, t, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
        	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
        	continuum_defln, hgtwnm, cenwnm, widwnm, fwnm, \
        	continuum_deflnyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x, \
			zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
			continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		nrgauss_cnm = len(taucnm1)
		nrg_cnm     = nrgauss_cnm

		x1 = x
		t1 = t
		# End 1665 BD1 #

		# 1667 BDBDBDDB 2 ============= #
		x = vlsr2
		t = ab_avg2

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x,zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
												continuum_defln, hgtwnm, cenwnm, widwnm, fwnm)

		tfita, sigma, \
		zrocnm2, taucnm2, cencnm2, widcnm2, tspincnm2, \
		sigzrocnm2, sigtaucnm2, sigcencnm2, sigwidcnm2, sigtspincnm2, \
		continuum_defln2, hgtwnm2, cenwnm2, widwnm2, fwnm2, \
		sigcontinuum_defln2, sighgtwnm2, sigcenwnm2, sigwidwnm2, sigfwnm2, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, x, t, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
        	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
        	continuum_defln, hgtwnm, cenwnm, widwnm, fwnm, \
        	continuum_deflnyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x, \
			zrocnm2, taucnm2, cencnm2, widcnm2, tspincnm2, ordercnm, \
			continuum_defln2, hgtwnm2, cenwnm2, widwnm2, fwnm2)

		nrgauss_cnm = len(taucnm2)
		nrg_cnm     = nrgauss_cnm

		x2 = x
		t2 = t

		# Tst6 #
		vplotrange= [velplotmin, velplotmax]
		te1= em_avg1
		te2= em_avg2

		tb_cont_1665, tb_wnm_tot, tb_cnm_tot, tb_tot_1665, exp_tau_sum = tb_exp(x1, \
			zrocnm1, taucnm1, cencnm1, \
			widcnm1, tspincnm1, ordercnm,\
			continuum_defln1, hgtwnm, cenwnm, widwnm, fwnm)

		tb_cont_1667, tb_wnm_tot, tb_cnm_tot, tb_tot_1667, exp_tau_sum = tb_exp(x2, \
			zrocnm1, taucnm1, cencnm1, \
			widcnm1, tspincnm1, ordercnm,\
			continuum_defln1, hgtwnm, cenwnm, widwnm, fwnm)

		# END - Tst6 #
		# ===============================#

		# ======== Do emission spectrum =================                                         
		# ;first do board 1
		x        = x1
		te       = te1

		tspincnm   = [4.]*npeak
		tspincnmyn = [1]*npeak
		tbaseline  = 89.05

		# Tst 7 #
		vrange     = [966, 1023]
		continuum_em = tbg

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

		tee = te - tbaseline + tbg
		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x, \
				zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
				continuum_em, hgtwnm, cenwnm, widwnm, fwnm)

		plt.plot(x,tee)
		plt.xlim(velplotmin, velplotmax)
		plt.plot(x,tb_cnm_tot+tbg)
		plt.plot(x, tb_tot)

		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_em1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, x, tee, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
	    	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
	    	continuum_em, hgtwnm, cenwnm, widwnm, fwnm, \
	    	continuum_emyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x, \
				zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
				continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		plt.plot(x, tb_tot)
		plt.plot(x,tee-tb_tot)

		nrgauss_wnm = len(hgtwnm1)
		nrg_wnm     = nrgauss_wnm
		# End  - Tst 7 #
		

		# Then do board 2 -- 1667
		x        = x2
		te       = te2

		tspincnm   = [4.]*npeak
		tspincnmyn = [1]*npeak
		tbaseline  = 90.55

		# Tst 7 #
		vrange       = [966, 1023]
		continuum_em = tbg

		zrocnm = zrocnm2
		taucnm = taucnm2
		cencnm = cencnm2
		widcnm = widcnm2
		nrcnm  = len(cencnm2)

		zrocnmyn       = 0
		continuum_emyn = 0
		taucnmyn       = [0]*nrcnm
		cencnmyn       = [0]*nrcnm
		widcnmyn       = [0]*nrcnm
		ordercnm       = list(range(nrcnm))

		tee = te - tbaseline + tbg
		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x, \
				zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
				continuum_em, hgtwnm, cenwnm, widwnm, fwnm)

		plt.plot(x,tee)
		plt.xlim(velplotmin, velplotmax)
		plt.plot(x,tb_cnm_tot+tbg)
		plt.plot(x, tb_tot)

		tfita, sigma, \
		zrocnm2, taucnm2, cencnm2, widcnm2, tspincnm2, \
		sigzrocnm2, sigtaucnm2, sigcencnm2, sigwidcnm2, sigtspincnm2, \
		continuum_em2, hgtwnm2, cenwnm2, widwnm2, fwnm2, \
		sigcontinuum_em2, sighgtwnm2, sigcenwnm2, sigwidwnm2, sigfwnm2, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, x, tee, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
	    	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
	    	continuum_em, hgtwnm, cenwnm, widwnm, fwnm, \
	    	continuum_emyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(x, \
				zrocnm2, taucnm2, cencnm2, widcnm2, tspincnm2, ordercnm, \
				continuum_em2, hgtwnm2, cenwnm2, widwnm2, fwnm2)

		plt.plot(x, tb_tot)
		plt.plot(x, tee-tb_tot)

		nrgauss_wnm = len(hgtwnm1)
		nrg_wnm     = nrgauss_wnm
		# End  - Tst 7 #

		print 'bd2...problem: ', problem

		print '####### 1665 ######'
		print '1) 1665 Tex & Error:'
		print '    ', tspincnm1
		print '    ', sigtspincnm1
		print '2) 1665 Tau & Error:'
		print '    ', taucnm1
		print '    ', sigtaucnm1
		print '3) 1665 V0 & Error:'
		print '    ', cencnm1
		print '    ', sigcencnm1
		print '4) 1665 Width & Error:'
		print '    ', widcnm1
		print '    ', sigwidcnm1
		print '5) 1665 hgtwnm & Error:'
		print '    ', hgtwnm1
		print '    ', sighgtwnm1
		print '6) 1665 cenwnm & Error:'
		print '    ', cenwnm1
		print '    ', sigcenwnm1
		print '7) 1665 widwnm & Error:'
		print '    ', widwnm1
		print '    ', sigwidwnm1
		print '7) 1665 fwnm & Error:'
		print '    ', fwnm1
		print '    ', sigfwnm1

		print ''
		print '####### 1667 ######'
		print '1) 1667 Tex & Error:'
		print '    ', tspincnm2
		print '    ', sigtspincnm2
		print '2) 1667 Tau & Error:'
		print '    ', taucnm2
		print '    ', sigtaucnm2
		print '3) 1667 V0 & Error:'
		print '    ', cencnm2
		print '    ', sigcencnm2
		print '4) 1667 Width & Error:'
		print '    ', widcnm2
		print '    ', sigwidcnm2
		print '5) 1667 hgtwnm & Error:'
		print '    ', hgtwnm2
		print '    ', sighgtwnm2
		print '6) 1667 cenwnm & Error:'
		print '    ', cenwnm2
		print '    ', sigcenwnm2
		print '7) 1667 widwnm & Error:'
		print '    ', widwnm2
		print '    ', sigwidwnm2
		print '7) 1667 fwnm & Error:'
		print '    ', fwnm2
		print '    ', sigfwnm2


        # plt.plot(xd,tb_tot, 'r-')
        # plt.xlim(xmin,xmax)
        # plt.grid()
        # plt.show()

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408
cal_tex(data, inf408)

sys.exit()