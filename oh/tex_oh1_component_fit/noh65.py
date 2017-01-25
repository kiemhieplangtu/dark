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

	em_avg1 = 0.5*correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = 0.5*correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = 0.5*correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = 0.5*correct_ctrl_chnl(data.la.i_abs_med_bd1)

	em_avg2 = 0.5*correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = 0.5*correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = 0.5*correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = 0.5*correct_ctrl_chnl(data.la.i_abs_med_bd2)

	# em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	# em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	# ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	# ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	# em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	# em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	# ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	# ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	return n,ell[n],bee[n],oh_f1[n],vlsr1[n],oh_f2[n],vlsr2[n],em_avg1[n],ab_avg1[n],em_avg2[n],ab_avg2[n]

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
def peak_guessp(sc,fname='../data/gauss_1665_peaks.txt'):
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

## Read Ningyu Tex##
 #
 # params string fname File-name
 # params string sc    Source

 # return list guessp Guess-Parameters 
 #        list base_range Range of Baselines
 #        float baseline Baselines after checking chi2
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_carl_tex(sc,fname='../sub_data/carl_tex.txt'):
	ret = {}

	cols = ['ng1','ng2', 'ng', 'src','l','b','zro1','tau1', 'cen1', 'wid1', 'tex1', 'tex1er', 'tex2', 'tex2er']
	fmt  = ['i',  'i',   'i',  's',  'f','f', 'f',   'f',   'f',     'f',    'f',   'f',     'f',      'f'  ]
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read()
	src  = info['src']
	tau  = info['tau1']
	tex  = info['tex1']
	for i in range(len(src)):
		if info['src'][i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['tau']        = [tau[i]]
			ret[src[i]]['tex65']      = [tex[i]]
		else:
			ret[src[i]]['tau']        += [tau[i]]
			ret[src[i]]['tex65']      += [tex[i]]

	return ret[sc]['tau'], ret[sc]['tex65']


## Compute Tex for 1665 line ##
 #
 # params dict data Data
 # params dict inf408  Info about the Tb_background at 408MHz
 #
 # return fit-results
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tex(data, inf408, bd=1):
	fit      = gfit()
	bg408    = read_tbg408_healpy('../result/bg408_to_compare.txt')
 	src_list = list(data.la.srcname)

 	for src in src_list:
	 	## Basic infor of the source
	 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src,src_list)

	 	## Vel range to cal. Baseline-background
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2 = read_bins_to_cal_bg(n)

	 	## Vel range to cal. Baseline-background
	 	xmin, xmax = vel_range(n)

	 	# if(src != '3C105'):
	 	# 	continue

	 	if(xmin==0 and xmax==0):
	 		continue

	 	## Absorption and Emission data ##
	 	v   = vlsr1
		Ta  = ab_avg1
		Te  = em_avg1
		cst = 3.99757843817
		frq = 1665.402
		pfl = '../data/gauss_1665_peaks.txt'
		if(bd == 2):
			v   = vlsr2
			Ta  = ab_avg2
			Te  = em_avg2
			cst = 2.21841824609
			frq = 1667.359
			pfl = '../data/gauss_1667_peaks.txt'

		## Ningyu Tex
		ntau, ntex = read_carl_tex(src,fname='../sub_data/carl_tex.txt')

		## Background ##
		# tbg1665    = 2.8+get_tb_408(ell,bee,inf408.tb_408)*(408./frq)**2.8 # Tbg from 408MHz
		# tbg1665    = 2.8+bg408[src]*(408./frq)**2.8 # Tbg from 408MHz

		## BACKGROUNDS and THEIR UNCERTAINTIES ##
		# tc1665,tc_er        = baseline_from_linear_fit(v, Ta, avmin1, avmax1, avmin2, avmax2,fit=False)
		# bg_off1665,bgoff_er = baseline_from_linear_fit(v, Te, evmin1, evmax1, evmin2, evmax2,fit=False)
		tau, v0, wid, base_range, v1,v2, tc1665, bg_off1665, tbg1665, base = peak_guessp(src,pfl)
		tc1665     = tc1665/2.
		bg_off1665 = bg_off1665/2.
		base       = base/2.
		trx        = bg_off1665 - tbg1665

		tbg = tbg1665
		continuum_em = tbg

		# Get Index of the Velocity-range #
		xmax_id  = get_vel_index(v, xmin)   # xmax_index
		xmin_id  = get_vel_index(v, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id          # Total number of bins
		vrange   = [xmin_id, xmax_id]
		vrange   = [v1, v2]
		dv       = (xmax-xmin)/num_chnl

		# Linear fit the baselines of spectra #
		# Vrange infor to calculate baseline - Linear fit #
		# print '***********************'
		# print n, ' ', data.la[n].srcname, data.la[n].ell, data.la[n].bee
		# print 'Tc: ', tc1665
		# print 'Emission baseline: ',bg_off1665
		# print 'Tbg_from_408: ',tbg1665	
		# print '-------- ***** ---------'

		npeak   = len(tau)
		zrolvyn = 1
		tauyn   = [1]*npeak
		v0yn    = [1]*npeak
		widyn   = [1]*npeak

		cont  = tc1665 #59.15
		zrolv = 0.
		tex   = [1.e-6]*npeak

		zrolvyn = 0
		tauyn   = [1]*npeak
		v0yn    = [1]*npeak
		widyn   = [1]*npeak
		texyn   = [0]*npeak
		contyn  = 1


		##  1665 ###
		tb_tot = fit.tb_exp(v,zrolv, tau, v0, wid, tex, cont)

		# plt.plot(v, Ta)
		# plt.plot(v, tb_tot, 'g-')
		# plt.title(data.la[n].srcname+' - 1665 Initial Fit')
		# plt.xlim(-20., 5.)
		# plt.grid()
		# plt.show()
		
		# Fit Absoprption line
		Tfita, sigma, \
		zrolvfit, taufit, v0fit, widfit, texfit, \
		zrolv_er, tau_er, v0_er, wid_er, tex_er, \
		contfit,\
		cont_er,\
		cov, nloop, nloopmax = fit.fit( v, Ta, vrange, \
			zrolv, tau, v0, wid, tex, \
	    	zrolvyn, tauyn, v0yn, widyn, texyn, \
	    	cont,contyn)

		tb_tot = fit.tb_exp(v, \
			zrolvfit, taufit, v0fit, widfit, texfit,\
			contfit)

		# plt.plot(v,Ta)
		# plt.plot(v, tb_tot, 'r-')
		# plt.title(data.la[n].srcname+' - 1665 Absortion line Fit')
		# plt.xlim(-20., 5.)
		# plt.grid()
		# plt.show()

		# ===============================#

		# ======== Emission Line =================
		tex   = [3.]*npeak
		texyn = [1]*npeak

		cont_em    = tbg
		tbaseline  = round(bg_off1665,3)
		Tefit      = Te - tbaseline + tbg

		zrolv = zrolvfit
		tau   = taufit
		v0    = v0fit
		wid   = widfit
		nrcnm = len(taufit)

		zrolvyn  = 0
		contemyn = 0
		tauyn    = [0]*npeak
		v0yn     = [0]*npeak
		widyn    = [0]*npeak

		tb_tot = fit.tb_exp(v, \
				zrolv, tau, v0, wid, tex, cont_em)

		# plt.plot(v,Tefit)
		
		# plt.plot(v,tb_tot,'g-')
		# plt.title(data.la[n].srcname+' - 1665 Emission First glance')
		# plt.xlim(-20., 5.)
		# plt.grid()
		# plt.show()

		# # Fit with different Baselines to get the best fit ##
		# # Uncomment to see the Error vs Baseline ##
		# bsl = []
		# er  = []
		# for x in np.arange(70.88,70.93,0.01):
		# 	tfita, sigma, \
		# 	zrolvfite1, taufite1, v0fite1, widfite1, texfite1, \
		# 	zrolv_ere1, tau_ere1, v0_ere1, wid_ere1, tex_ere1, \
		# 	cont_eme1,\
		# 	cont_ere1,\
		# 	cove1, nloope1, nloopmax = fit.fit( v, Te - x + tbg, vrange, \
		# 		zrolv, tau, v0, wid, tex, \
		#     	zrolvyn, tauyn, v0yn, widyn, texyn, \
		#     	cont_em,\
		#     	contemyn)
		# 	bsl.append(x)
		# 	er.append(sigma)

		# tb_tot = fit.tb_exp(v, \
		# 		zrolvfite1, taufite1, v0fite1, widfite1, texfite1, cont_eme1)

		# plt.plot(v,Tefit)
		# plt.plot(v, tb_tot,'r-')
		# plt.title(data.la[n].srcname+' - 1665 Emission fit')
		# plt.xlim(-20., 5.)
		# plt.show()

		# plt.plot(bsl,er)
		# plt.title(data.la[n].srcname+' - Error vs baseline')
		# plt.show()

		## After check Chi2 vs baseline, Select the baseline with the best fit ##
		tbaseline  = base #70.90
		Tefit      = Te - tbaseline + tbg

		Tfite, sigma, \
		zrolvfite, taufite, v0fite, widfite, texfite, \
		zrolv_ere, tau_ere, v0_ere, wid_ere, tex_ere, \
		cont_eme,\
		cont_ere,\
		cove, nloope, nloopmax = fit.fit(v, Tefit, vrange, \
			zrolv, tau, v0, wid, tex, \
	    	zrolvyn, tauyn, v0yn, widyn, texyn, \
	    	cont_em,\
	    	contemyn)

		# plt.plot(v,Tefit)
		# plt.plot(v, tb_tot,'r-')
		# plt.title(data.la[n].srcname+' - 1665 Emission Final fit')
		# plt.xlim(-20., 5.)
		# plt.grid()
		# plt.show()

		texfit = texfite
		tex_er = tex_ere

	    ## Calculate N(OH) ##
		popt = []
		for i in range(len(tau)):
			popt += [tau[i], v0[i], wid[i],texfit[i]]

		stau_fit = quad(tau_func, xmin,xmax, args=tuple(popt))
		# print '0) Tex[i]*integral(tau)'
		# print '    ', stau_fit

		noh_fit = 3.99757843817*stau_fit[0]  # x10^14

		# print '1) 1665 Tex & Error:'
		# print '    ', texfit
		# print '    ', tex_ere
		# print '2) Baseline final fit:'
		# print '    ', tbaseline
		# print '3) 1665 Tau & Error:'
		# print '    ', taufit
		# print '    ', tau_er
		# print '3) 1665 V0 & Error:'
		# print '    ', v0
		# print '    ', v0_er
		# print '3) 1665 Width & Error:'
		# print '    ', widfit
		# print '    ', wid_er
		# print '4) N(OH):'
		# print '    ', noh_fit, 'x10^14'
		# print '5) N loops for Emission line'
		# print '    ', nloope

		for k in range(npeak):
			tauk    = round(taufit[k],8)
			ntauk   = round(ntau[k],8)
			tauerk  = round(tau_er[k],8)
			v0k     = round(v0[k],8)
			v0erk   = round(v0_er[k],8)
			widk    = round(widfit[k],8)
			widerk  = round(wid_er[k],8)
			texk    = round(texfit[k],6)
			ntexk   = ntex[k]
			tex_erk = round(tex_er[k],6)
			print '{}  {}\t{}\t{:04.8f}\t{:04.6f}\t{:04.6f}\t{:04.6f}\t{:04.6f}\t{}\t{:04.6f}'\
						.format(n, src, str(tauk)+'/'+str(ntauk), tauerk, v0k, v0erk, widk, widerk, str(texk)+'/'+ str(ntexk), tex_erk )

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408
cal_tex(data, inf408)

sys.exit()