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
def gcurv(xdata, zro1, hgt1, cen1, widfit):
	#DETERMINE NR OF GAUSSIANS...
	ngaussians = 1 # cho vui

	tfit = 0.*xdata + zro1
	if (widfit > 0.):
		tfit = tfit + hgt1*np.exp(- ( (xdata-cen1)/(0.6005612*widfit))**2)

	return tfit

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

	src  = '3C18'
	n    = 18
	
	xmin = vmin[n]
	xmax = vmax[n]

	tbg = 3.44514
	continuum_em = tbg

	# VLSR #
	x = vlsr1[n]
	velplotmin = xmin
	velplotmax = xmax
	# Get Index of the Velocity-range #
	xmax_id  = get_vel_index(x, xmin)   # xmax_index
	xmin_id  = get_vel_index(x, xmax)   # xmin_index
	num_chnl = xmax_id-xmin_id          # Total number of bins
	vrange   = [xmin_id, xmax_id]
	vrange   = [879,1460]

	print n, '   ', data.la[n].srcname, data.la[n].ell, data.la[n].bee
	zrolvyn = 1
	tauyn   = [1,1]
	v0yn    = [1,1]
	widyn   = [1,1]

	cont  = 59.15
	zrolv = 0.
	tau   = [0.1/cont, 0.6/cont]
	v0    = [-9.9, -7.6]
	wid   = [1.3, 0.9]
	tex   = [1.e-6, 1.e-6]

	zrolvyn = 0
	tauyn   = [1, 1]
	v0yn    = [1, 1]
	widyn   = [1,1]
	texyn   = [0, 0]
	contyn  = 1

	##  1665  bd1 ###
	##  1665  bd1 ###
	##  1665  bd1 ###
	xd = vlsr1[n]
	td = ab_avg1[n]

	xd_bd1 = xd
	td_bd1 = td

	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd,zrolv, tau, v0, wid, tex, cont)

	plt.plot(xd, td)
	plt.plot(xd, tb_cont, 'g-')
	plt.title('1665 Initial Fit')
	plt.xlim(-20., 5.)
	plt.show()
	
	# Fit Absoprption line
	tfita, sigma, \
	zrolvfit1, taufit1, v0fit1, widfit1, texfit1, \
	zrolv_er1, tau_er1, v0_er1, wid_er1, tex_er1, \
	contfit1,\
	cont_er1,\
	cov1, nloop1, \
	tb_cont1, tb_cnm_tot1, \
	exp_tau_sum1, nloopmax = fit.fit( xd, td, vrange, \
		zrolv, tau, v0, wid, tex, \
    	zrolvyn, tauyn, v0yn, widyn, texyn, \
    	cont,contyn)

	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd, \
		zrolvfit1, taufit1, v0fit1, widfit1, texfit1,\
		contfit1)

	plt.plot(xd,td)
	plt.plot(xd, tb_cont, 'r-')
	plt.title('1665 Absortion line Fit')
	plt.xlim(-20., 5.)
	plt.show()
	# End 1665 BD1 #

	# 1667 BDBDBDDB 2 ============= #
	xd = vlsr2[n]
	td = ab_avg2[n]

	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd,zrolv, tau, v0, wid, tex, cont)

	tfita, sigma, \
	zrolvfit2, taufit2, v0fit2, widfit2, texfit2, \
	zrolv_er2, tau_er2, v0_er2, wid_er2, tex_er2, \
	contfit2, cont_er2,\
	cov2, nloop2, \
	tb_cont2, tb_cnm_tot2, \
	exp_tau_sum2, nloopmax = fit.fit( xd, td, vrange, \
		zrolv, tau, v0, wid, tex, \
    	zrolvyn, tauyn, v0yn, widyn, texyn, \
    	cont,\
    	contyn)

	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd, \
		zrolvfit2, taufit2, v0fit2, widfit2, texfit2, contfit2)

	xd_bd2 = xd
	td_bd2 = td

	# Tst6 #
	tde_bd1 = em_avg1[n]
	tde_bd2 = em_avg2[n]

	# END - Tst6 #
	# ===============================#

	# ======== Emission Spectra =================                                         
	# ;first do 1665  
	xd    = xd_bd1
	tde   = tde_bd1

	tex   = [4.,7.]
	texyn = [1,1]
	tbaseline  = 70.90

	# Tst 7 #
	cont_em = tbg

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

	tdee = tde - tbaseline + tbg
	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd, \
			zrolv, tau, v0, wid, tex, cont_em)

	plt.plot(xd,tdee)
	
	plt.plot(xd,tb_cnm_tot+tbg,'r-')
	plt.title('1665 Emission First glance')
	plt.xlim(-20., 5.)
	plt.show()

	tfita, sigma, \
	zrolvfite1, taufite1, v0fite1, widfite1, texfite1, \
	zrolv_ere1, tau_ere1, v0_ere1, wid_ere1, tex_ere1, \
	cont_eme1,\
	cont_ere1,\
	cove1, nloope1, \
	tb_conte1, tb_cnm_tote1, \
	exp_tau_sume1, nloopmax = fit.fit( xd, tdee, vrange, \
		zrolv, tau, v0, wid, tex, \
    	zrolvyn, tauyn, v0yn, widyn, texyn, \
    	cont_em,\
    	contemyn)

	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd, \
			zrolvfite1, taufite1, v0fite1, widfite1, texfite1, cont_eme1)

	plt.plot(xd,tdee)
	plt.plot(xd, tb_tot,'r-')
	plt.title('1665 Emission fit')
	plt.xlim(-20., 5.)
	plt.show()

	texfit1 = texfite1 ## DAU DAY
	tex_er1 = tex_ere1 ## DAU DAY

	#### OH 1667 ###
	xd       = xd_bd2
	tde      = tde_bd2

	tex   = [4.,1.]
	texyn = [1,0]
	tbaseline  = 71.40

	# Tst 7 #
	cont_em = tbg

	zrolv = zrolvfit2
	tau   = taufit2
	v0    = v0fit2
	wid   = widfit2
	nrcnm = len(taufit2)

	zrolvyn   = 0
	cont_emyn = 0
	tauyn     = [0]*nrcnm
	v0yn      = [0]*nrcnm
	widyn     = [0]*nrcnm

	tdee = tde - tbaseline + tbg
	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd, \
			zrolv, tau, v0, wid, tex, \
			continuum_em)

	plt.plot(xd, tdee)
	plt.plot(xd, tb_tot)
	plt.title('1667 Emission 1st glance')
	plt.xlim(-20., 5.)
	plt.show()

	tfita, sigma, \
	zrolvfite2, taufite2, v0fite2, widfite2, texfite2, \
	zrolv_ere2, tau_ere2, v0_ere2, wid_ere2, tex_ere2, \
	cont_eme2,\
	cont_emer2,\
	cove2, nloope2, \
	tb_conte2, tb_cnm_tote2, \
	exp_tau_sume2, nloopmax = fit.fit( xd, tdee, vrange, \
		zrolv, tau, v0, wid, tex, \
    	zrolvyn, tauyn, v0yn, widyn, texyn, \
    	cont_em,\
    	cont_emyn)

	tb_cont, tb_cnm_tot, tb_tot, exp_tau_sum = fit.tb_exp(xd, \
			zrolvfite2, taufite2, v0fite2, widfite2, texfite2, \
			cont_eme2)

	plt.plot(xd, tdee)
	plt.plot(xd, tb_tot)
	plt.title('1667 Emission fit')
	plt.xlim(-20., 5.)
	plt.show()

	texfit2 = texfite2 ## DAU DAY
	tex_er2 = tex_er2 ## DAU DAY
	# End  - Tst 7 #

	print '1665 Ex temps:'
	print texfit1
	print tex_ere1
	print nloope1

	print '1667 Ex temps:'
	print texfit2
	print tex_er2
	print nloope2

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408
cal_tex(data, inf408)

sys.exit()