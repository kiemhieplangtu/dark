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

## Compute the Tex spectrum from etau ##
 #
 # params list etau_fit Exp(-tau)
 # params list velo     List of Velocity
 # params list tb_off   List of Off-source Temperature
 # params float tbg     CMB + Radio noise
 # params float trx     Receiver Temperature
 # params list popt     Fit-result parameters 
 #
 # return List tex      Excitation temperature
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tex_from_etau(etau_fit,velo,tb_off,tbg,trx,popt):
	tex  = []

	# Find the vel-range of the peak #
	peak = cal_vel_range_fwhm(popt)

	tex_peak  = []
	tex_count = []
	for k in range(0,len(peak)/2):
		tex_peak.append(0.)
		tex_count.append(0)

	# Compute the Excitation Temperature #
	s     = 0.
	count = 0
	for i in range(0, len(etau_fit)):
		if (etau_fit[i] != 1.):
			ts = (tb_off[i]-trx-tbg*etau_fit[i])/(1.-etau_fit[i])
			tex.append(ts)
		else:
			tex.append(0.)

		for k in range(0,len(peak),2):
			vmin = peak[0+k]
			vmax = peak[1+k]

			if ((velo[i]>vmin) and (velo[i]<vmax)) :
				tex_peak[k/2]  = tex_peak[k/2] + tex[i]
				tex_count[k/2] = tex_count[k/2] + 1

				s     = s + tex[i]
				count = count + 1

	print 'Tex for Peaks'
	ss = 0.
	for k in range(0,len(peak)/2):
		tex_peak[k]  = tex_peak[k]/tex_count[k]
		ss           = ss + tex_peak[k]
		print tex_peak[k],tex_count[k]
	print 'TEX ==== '
	print s/count, ss, count
	print 'TEX ==== '
	s = s/count
	return tex,s

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

## Test ##
 #
 # params 
 # params 
 #
 # return void
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def tst7(look,xrange,velplotmin, velplotmax,amg_1666,\
	nrcnm,tde,tbaseline,tbg,xd,zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm,\
	continuum_em, hgtwnm, cenwnm, widwnm, fwnm,\
	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn,\
	continuum_emyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn):

	# Tst 7 #
	zrocnmyn       = 0
	continuum_emyn = 0
	taucnmyn       = [0]*nrcnm
	cencnmyn       = [0]*nrcnm
	widcnmyn       = [0]*nrcnm
	ordercnm       = [0]*nrcnm

	tdee = tde - tbaseline + tbg
	tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
			continuum_em, hgtwnm, cenwnm, widwnm, fwnm)

	plt.plot(xd,tdee)
	plt.xlim(velplotmin, velplotmax)
	plt.plot(xd,tb_cnm_tot+tbg)
	plt(xd, tb_tot)
	# plt.show()

	tfita, sigma, \
	zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
	sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
	continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
	sigcontinuum_em1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
	cov, problem, nloop, \
	tb_cont, tb_wnm_tot, tb_cnm_tot, \
	exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, xd, td, vrange, \
		zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
    	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
    	continuum_defln, hgtwnm, cenwnm, widwnm, fwnm, \
    	continuum_deflnyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)


	tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
			zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
			continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

	plt.plot(xd, tb_tot)
	plt.plot(xd,tdee-tb_tot)

	nrgauss_wnm = len(hgtwnm1)
	nrg_wnm     = nrgauss_wnm
	emg_1666    = {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
	      'tbaseline':0.0, \
	      'continuum':0.0, 'hgtwnm':0.0, 'cenwnm':0.0, 'widwnm':0.0, 'fwnm':0.0, \
	      'sigcontinuum':0.0, 'sighgtwnm':0.0, 'sigcenwnm':0.0, 'sigwidwnm':0.0, 'sigfwnm':0.0,\
	      'continuum_em':0.0}

	dassign(emg_1666,'sname', src, nrgauss_wnm)
	dassign(emg_1666,'ell', ell[n], nrgauss_wnm)
	dassign(emg_1666,'bee', bee[n], nrgauss_wnm)
	dassign(emg_1666,'gaussnr', list(range(nrg_wnm)))
	dassign(amg_1665,'nrgauss', nrg_wnm, nrgauss_wnm)

	dassign(emg_1666,'tbaseline', tbaseline, nrgauss_wnm)
	dassign(emg_1666,'continuum_em', continuum_em, nrgauss_wnm)

	dassign(emg_1666,'hgtwnm', hgtwnm1)
	dassign(emg_1666,'cenwnm', cenwnm1)
	dassign(emg_1666,'widwnm', widwnm1)
	dassign(emg_1666,'fwnm', fwnm1)

	dassign(emg_1666,'sighgtwnm', sighgtwnm1)
	dassign(emg_1666,'sigcenwnm', sigcenwnm1)
	dassign(emg_1666,'sigwidwnm', sigwidwnm1)
	dassign(emg_1666,'sigfwnm', sigfwnm1)

	amg_1666['tspincnm']    = tspincnm1 ## DAU DAY
	amg_1666['sigtspincnm'] = sigtspincnm1 ## DAU DAY

	# End  - Tst 7 #

	return ''

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
 	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee


	oh_f = data.la.cfr_bd1
	vlsr = data.la.vlsr_bd1

	em_avg = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	# Gaussian peaks' info - estimate values of Amp, tau0 and Width #
	peak = peak_info('data/gauss_1665_peaks.txt')

	# Vrange infor #
	fname    = 'data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']

	# 26 src with no CO #
	fname    = 'data/26src_no_co.txt'
	cols     = ['idx','src']
	fmt      = ['i','s']

	src_no_co = restore(fname, 2, cols, fmt)
	s26info   = src_no_co.read()
	s26src    = s26info['src']

	amgoh= {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
      'zrocnm':0, 'taucnm':0.0, 'cencnm':0.0, 'widcnm':0.0, 'tspincnm':0.0, \
      'sigzrocnm':0, 'sigtaucnm':0.0, 'sigcencnm':0.0, 'sigwidcnm':0.0, 'sigtspincnm':0.0, \
      'continuum_defln':0.0, 'sigcontinuum_defln':0.0, 'continuum_em':0.0}

	emgoh= {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
	      'tbaseline':0.0, \
	      'continuum':0.0, 'hgtwnm':0.0, 'cenwnm':0.0, 'widwnm':0.0, 'fwnm':0.0, \
	      'sigcontinuum':0.0, 'sighgtwnm':0.0, 'sigcenwnm':0.0, 'sigwidwnm':0.0, 'sigfwnm':0.0,\
	      'continuum_em':0.0}

	for src in src_list:
		if(src != '3C123'):
			continue

		n    = src_list.index(src)
		
		data.la[n].i_abs_avg_bd0[1023]= (data.la[n].i_abs_avg_bd0[1022]+ data.la[n].i_abs_avg_bd0[1024])/2.
		data.la[n].i_abs_avg_bd1[1023]= (data.la[n].i_abs_avg_bd1[1022]+ data.la[n].i_abs_avg_bd1[1024])/2.
		data.la[n].i_abs_avg_bd2[1023]= (data.la[n].i_abs_avg_bd2[1022]+ data.la[n].i_abs_avg_bd2[1024])/2.

		data.la[n].i_em_avg_bd0[1023] = (data.la[n].i_em_avg_bd0[1022]+ data.la[n].i_em_avg_bd0[1024])/2.
		data.la[n].i_em_avg_bd1[1023] = (data.la[n].i_em_avg_bd1[1022]+ data.la[n].i_em_avg_bd1[1024])/2.
		data.la[n].i_em_avg_bd2[1023] = (data.la[n].i_em_avg_bd2[1022]+ data.la[n].i_em_avg_bd2[1024])/2.

		
		xmin = vmin[n]
		xmax = vmax[n]

		if (xmin == 0. and xmax == 0.):
			continue

		tbg=50./2
		continuum_em= tbg

		# VLSR #
		x = vlsr[n]
		velplotmin = xmin
		velplotmax = xmax
		# Get Index of the Velocity-range #
		xmax_id  = get_vel_index(x, xmin)   # xmax_index
		xmin_id  = get_vel_index(x, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id          # Total number of bins
		vrange   = [xmin_id, xmax_id]
		vrange   = [966,1023]

		print n, '   ', data.la[n].srcname, data.la[n].ell, data.la[n].bee

		look   = -1
		zro0yn = 1
		hgt0yn = [1,1,1]
		cen0yn = [1,1,1]
		wid0yn = [1,1,1]

		continuum_defln = 532.
		continuumyn = 1

		zrocnm   = 0.
		taucnm   = [30./continuum_defln, 30./continuum_defln, 30./continuum_defln]
		cencnm   = [3.2, 4.4, 5.3]
		widcnm   = [0.6,0.6,0.6]
		tspincnm = [1.e-6, 1.e-6,1.e-6]
		ordercnm = [0,1,2]

		zrocnmyn   = 0
		taucnmyn   = [1, 1,1]
		cencnmyn   = [1, 1,1]
		widcnmyn   = [1,1,1]
		tspincnmyn = [0, 0,0]

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
		xd = data.la[n].vlsr_bd1
		td = data.la[n].i_abs_avg_bd1

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd,zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
												continuum_defln, hgtwnm, cenwnm, widwnm, fwnm)

		plt.plot(xd,td)
		plt.plot(xd, tb_cont, 'g-')

		fit  = gfit()
		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_defln1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, xd, td, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
        	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
        	continuum_defln, hgtwnm, cenwnm, widwnm, fwnm, \
        	continuum_deflnyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
			zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
			continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		nrgauss_cnm = len(taucnm1)
		nrg_cnm     = nrgauss_cnm
		amg_1665    = {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
      		'zrocnm':0, 'taucnm':0.0, 'cencnm':0.0, 'widcnm':0.0, 'tspincnm':0.0, \
      		'sigzrocnm':0, 'sigtaucnm':0.0, 'sigcencnm':0.0, 'sigwidcnm':0.0, 'sigtspincnm':0.0, \
      		'continuum_defln':0.0, 'sigcontinuum_defln':0.0, 'continuum_em':0.0}
		dassign(amg_1665,'sname', src, nrgauss_cnm)
		dassign(amg_1665,'ell', ell[n], nrgauss_cnm)
		dassign(amg_1665,'bee', bee[n], nrgauss_cnm)
		dassign(amg_1665,'gaussnr', list(range(nrg_cnm)))
		dassign(amg_1665,'nrgauss', nrg_cnm, nrgauss_cnm)

		dassign(amg_1665,'zrocnm', zrocnm1, nrgauss_cnm)
		dassign(amg_1665,'taucnm', taucnm1)
		dassign(amg_1665,'cencnm', cencnm1)
		dassign(amg_1665,'widcnm', widcnm1)
		dassign(amg_1665,'tspincnm', tspincnm1)

		dassign(amg_1665,'sigzrocnm', sigzrocnm1, nrgauss_cnm)
		dassign(amg_1665,'sigtaucnm', sigtaucnm1.tolist())
		dassign(amg_1665,'sigcencnm', sigcencnm1.tolist())
		dassign(amg_1665,'sigwidcnm', sigwidcnm1)
		dassign(amg_1665,'sigtspincnm', sigtspincnm1.tolist())

		dassign(amg_1665,'continuum_defln', continuum_defln1, nrgauss_cnm)
		dassign(amg_1665,'sigcontinuum_defln', sigcontinuum_defln1, nrgauss_cnm)
		dassign(amg_1665,'continuum_em', continuum_em, nrgauss_cnm)

		xd_bd1 = xd
		td_bd1 = td
		# End 1665 BD1 #

		# 1667 BDBDBDDB 2 ============= #
		xd = data.la[n].vlsr_bd2
		td = data.la[n].i_abs_avg_bd2

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd,zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
												continuum_defln, hgtwnm, cenwnm, widwnm, fwnm)

		fit  = gfit()
		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_defln1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, xd, td, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
        	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
        	continuum_defln, hgtwnm, cenwnm, widwnm, fwnm, \
        	continuum_deflnyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
			zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
			continuum_defln1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		nrgauss_cnm = len(taucnm1)
		nrg_cnm     = nrgauss_cnm
		amg_1667    = {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
      		'zrocnm':0, 'taucnm':0.0, 'cencnm':0.0, 'widcnm':0.0, 'tspincnm':0.0, \
      		'sigzrocnm':0, 'sigtaucnm':0.0, 'sigcencnm':0.0, 'sigwidcnm':0.0, 'sigtspincnm':0.0, \
      		'continuum_defln':0.0, 'sigcontinuum_defln':0.0, 'continuum_em':0.0}
		dassign(amg_1667,'sname', src, nrgauss_cnm)
		dassign(amg_1667,'ell', ell[n], nrgauss_cnm)
		dassign(amg_1667,'bee', bee[n], nrgauss_cnm)
		dassign(amg_1667,'gaussnr', list(range(nrg_cnm)))
		dassign(amg_1665,'nrgauss', nrg_cnm, nrgauss_cnm)

		dassign(amg_1667,'zrocnm', zrocnm1, nrgauss_cnm)
		dassign(amg_1667,'taucnm', taucnm1)
		dassign(amg_1667,'cencnm', cencnm1)
		dassign(amg_1667,'widcnm', widcnm1)
		dassign(amg_1667,'tspincnm', tspincnm1)

		dassign(amg_1667,'sigzrocnm', sigzrocnm1, nrgauss_cnm)
		dassign(amg_1667,'sigtaucnm', sigtaucnm1.tolist())
		dassign(amg_1667,'sigcencnm', sigcencnm1.tolist())
		dassign(amg_1667,'sigwidcnm', sigwidcnm1)
		dassign(amg_1667,'sigtspincnm', sigtspincnm1.tolist())

		dassign(amg_1667,'continuum_defln', continuum_defln1, nrgauss_cnm)
		dassign(amg_1667,'sigcontinuum_defln', sigcontinuum_defln1, nrgauss_cnm)
		dassign(amg_1667,'continuum_em', continuum_em, nrgauss_cnm)

		xd_bd2 = xd
		td_bd2 = td

		# Tst6 #
		vplotrange= [velplotmin, velplotmax]
		tde_bd1= data.la[n].i_em_avg_bd1
		tde_bd2= data.la[n].i_em_avg_bd2

		tb_cont_1665, tb_wnm_tot, tb_cnm_tot, tb_tot_1665, exp_tau_sum = tb_exp(xd_bd1, \
			amg_1665['zrocnm'][0], amg_1665['taucnm'], amg_1665['cencnm'], \
			amg_1665['widcnm'], amg_1665['tspincnm'], ordercnm,\
			amg_1665['continuum_defln'][0], hgtwnm, cenwnm, widwnm, fwnm)

		tb_cont_1667, tb_wnm_tot, tb_cnm_tot, tb_tot_1667, exp_tau_sum = tb_exp(xd_bd2, \
			amg_1667['zrocnm'][0], amg_1667['taucnm'], amg_1667['cencnm'], \
			amg_1667['widcnm'], amg_1667['tspincnm'], ordercnm,\
			amg_1667['continuum_defln'][0], hgtwnm, cenwnm, widwnm, fwnm)

		# END - Tst6 #
		# ===============================#

		# ======== do emission spectrum =================                                         
		# ;first do board 1   
		amg_1666 = amg_1665
		xd       = xd_bd1
		tde      = tde_bd1

		tspincnm   = [4.,4.,4.]
		tspincnmyn = [1,1,1]
		tbaseline  = 89.05

		# Tst 7 #
		vrange     = [966, 1023]
		continuum_em = tbg

		zrocnm = amg_1666['zrocnm'][0]
		taucnm = amg_1666['taucnm']
		cencnm = amg_1666['cencnm']
		widcnm = amg_1666['widcnm']
		nrcnm  = len(amg_1666['cencnm'])

		zrocnmyn       = 0
		continuum_emyn = 0
		taucnmyn       = [0]*nrcnm
		cencnmyn       = [0]*nrcnm
		widcnmyn       = [0]*nrcnm
		ordercnm       = list(range(nrcnm))

		tdee = tde - tbaseline + tbg
		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
				zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
				continuum_em, hgtwnm, cenwnm, widwnm, fwnm)

		plt.plot(xd,tdee)
		plt.xlim(velplotmin, velplotmax)
		plt.plot(xd,tb_cnm_tot+tbg)
		plt.plot(xd, tb_tot)

		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_em1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, xd, tdee, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
	    	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
	    	continuum_em, hgtwnm, cenwnm, widwnm, fwnm, \
	    	continuum_emyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
				zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
				continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		plt.plot(xd, tb_tot)
		plt.plot(xd,tdee-tb_tot)

		nrgauss_wnm = len(hgtwnm1)
		nrg_wnm     = nrgauss_wnm
		emg_1666    = {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
		      'tbaseline':0.0, \
		      'continuum':0.0, 'hgtwnm':0.0, 'cenwnm':0.0, 'widwnm':0.0, 'fwnm':0.0, \
		      'sigcontinuum':0.0, 'sighgtwnm':0.0, 'sigcenwnm':0.0, 'sigwidwnm':0.0, 'sigfwnm':0.0,\
		      'continuum_em':0.0}

		dassign(emg_1666,'sname', src, nrgauss_wnm)
		dassign(emg_1666,'ell', ell[n], nrgauss_wnm)
		dassign(emg_1666,'bee', bee[n], nrgauss_wnm)
		dassign(emg_1666,'gaussnr', list(range(nrg_wnm)))
		dassign(emg_1666,'nrgauss', nrg_wnm, nrgauss_wnm)

		dassign(emg_1666,'tbaseline', tbaseline, nrgauss_wnm)
		dassign(emg_1666,'continuum_em', continuum_em, nrgauss_wnm)

		dassign(emg_1666,'hgtwnm', hgtwnm1)
		dassign(emg_1666,'cenwnm', cenwnm1)
		dassign(emg_1666,'widwnm', widwnm1)
		dassign(emg_1666,'fwnm', fwnm1)

		dassign(emg_1666,'sighgtwnm', sighgtwnm1)
		dassign(emg_1666,'sigcenwnm', sigcenwnm1)
		dassign(emg_1666,'sigwidwnm', sigwidwnm1)
		dassign(emg_1666,'sigfwnm', sigfwnm1)

		amg_1666['tspincnm']    = tspincnm1 ## DAU DAY
		amg_1666['sigtspincnm'] = sigtspincnm1 ## DAU DAY
		# End  - Tst 7 #

		amg_1665 = amg_1666
		emg_1665 = emg_1666
		

		# ;Then do board 2 -- 1667   
		amg_1666 = amg_1667
		xd       = xd_bd2
		tde      = tde_bd2

		tspincnm   = [4.,4.,4.]
		tspincnmyn = [1,1,1]
		tbaseline  = 90.55

		# Tst 7 #
		vrange     = [966, 1023]
		continuum_em = tbg

		zrocnm = amg_1666['zrocnm'][0]
		taucnm = amg_1666['taucnm']
		cencnm = amg_1666['cencnm']
		widcnm = amg_1666['widcnm']
		nrcnm  = len(amg_1666['cencnm'])

		zrocnmyn       = 0
		continuum_emyn = 0
		taucnmyn       = [0]*nrcnm
		cencnmyn       = [0]*nrcnm
		widcnmyn       = [0]*nrcnm
		ordercnm       = list(range(nrcnm))

		tdee = tde - tbaseline + tbg
		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
				zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
				continuum_em, hgtwnm, cenwnm, widwnm, fwnm)

		plt.plot(xd,tdee)
		plt.xlim(velplotmin, velplotmax)
		plt.plot(xd,tb_cnm_tot+tbg)
		plt.plot(xd, tb_tot)

		tfita, sigma, \
		zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, \
		sigzrocnm1, sigtaucnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
		continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
		sigcontinuum_em1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
		cov, problem, nloop, \
		tb_cont, tb_wnm_tot, tb_cnm_tot, \
		exp_tau_sum, nloopmax, halfasseduse = fit.fit(look, xd, tdee, vrange, \
			zrocnm, taucnm, cencnm, widcnm, tspincnm, ordercnm, \
	    	zrocnmyn, taucnmyn, cencnmyn, widcnmyn, tspincnmyn, \
	    	continuum_em, hgtwnm, cenwnm, widwnm, fwnm, \
	    	continuum_emyn, hgtwnmyn, cenwnmyn, widwnmyn, fwnmyn)

		tb_cont, tb_wnm_tot, tb_cnm_tot, tb_tot, exp_tau_sum = tb_exp(xd, \
				zrocnm1, taucnm1, cencnm1, widcnm1, tspincnm1, ordercnm, \
				continuum_em1, hgtwnm1, cenwnm1, widwnm1, fwnm1)

		plt.plot(xd, tb_tot)
		plt.plot(xd,tdee-tb_tot)

		nrgauss_wnm = len(hgtwnm1)
		nrg_wnm     = nrgauss_wnm
		emg_1666    = {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
		      'tbaseline':0.0, \
		      'continuum':0.0, 'hgtwnm':0.0, 'cenwnm':0.0, 'widwnm':0.0, 'fwnm':0.0, \
		      'sigcontinuum':0.0, 'sighgtwnm':0.0, 'sigcenwnm':0.0, 'sigwidwnm':0.0, 'sigfwnm':0.0,\
		      'continuum_em':0.0}

		dassign(emg_1666,'sname', src, nrgauss_wnm)
		dassign(emg_1666,'ell', ell[n], nrgauss_wnm)
		dassign(emg_1666,'bee', bee[n], nrgauss_wnm)
		dassign(emg_1666,'gaussnr', list(range(nrg_wnm)))
		dassign(emg_1666,'nrgauss', nrg_wnm, nrgauss_wnm)

		dassign(emg_1666,'tbaseline', tbaseline, nrgauss_wnm)
		dassign(emg_1666,'continuum_em', continuum_em, nrgauss_wnm)

		dassign(emg_1666,'hgtwnm', hgtwnm1)
		dassign(emg_1666,'cenwnm', cenwnm1)
		dassign(emg_1666,'widwnm', widwnm1)
		dassign(emg_1666,'fwnm', fwnm1)

		dassign(emg_1666,'sighgtwnm', sighgtwnm1)
		dassign(emg_1666,'sigcenwnm', sigcenwnm1)
		dassign(emg_1666,'sigwidwnm', sigwidwnm1)
		dassign(emg_1666,'sigfwnm', sigfwnm1)

		amg_1666['tspincnm']    = tspincnm1 ## DAU DAY
		amg_1666['sigtspincnm'] = sigtspincnm1 ## DAU DAY
		# End  - Tst 7 #

		amg_1667 = amg_1666
		emg_1667 = emg_1666
		print 'bd2...problem: ', problem

		print '1665 spin temps:'
		print amg_1665['tspincnm']
		print amg_1665['sigtspincnm']
		print nloop

		print '1667 spin temps:'
		print amg_1667['tspincnm']
		print amg_1667['sigtspincnm']
		print nloop

		print
		print 'WNM 1665:'
		print emg_1665['hgtwnm']
		print emg_1665['sighgtwnm']
		print nloop

		print 'WNM 1667:'
		print emg_1667['hgtwnm']
		print emg_1667['sighgtwnm']
		print nloop






        # plt.plot(xd,tb_tot, 'r-')
        # plt.xlim(xmin,xmax)
        # plt.grid()
        # plt.show()

#============== MAIN ==============#
data   = readsav('data/makelines.sav') #data.la
inf408 = readsav('data/tb_408.sav') # l_cntr, b_cntr, tb_408
cal_tex(data, inf408)

sys.exit()