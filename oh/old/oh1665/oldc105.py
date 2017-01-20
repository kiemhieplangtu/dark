import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot
from gfit2          import gfit

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

## Test ##
 #
 # params 
 # params 
 #
 # return void
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def gauss_fit_step1(cfit,xd,tde,vrange,amg_oh,tbaseline,contoff,ts,tsyn):
	cont_em = contoff ## tbg=50./2

	tbg    = amg_oh['tbg'][0]
	tau    = amg_oh['tau']
	v0     = amg_oh['v0']
	wid    = amg_oh['wid']
	npeaks = len(amg_oh['v0'])

	tbgyn  = 0
	contyn = 0
	tauyn  = [0]*npeaks
	v0yn   = [0]*npeaks
	widyn  = [0]*npeaks
	order  = list(range(npeaks))

	tdee   = tde - tbaseline + contoff
	tb_tot = cfit.get_tb(xd, tbg, tau, v0, wid, ts, cont_em)

	tfit, error, \
	tbgfit, taufit, v0fit, widfit, \
	tbgerr, tauerr, v0err, widerr, \
	cont,conterr,\
	cov, problem, nloop, \
	nloopmax = cfit.fit(xd,tdee,vrange,cont_em,tbg,tau,v0,wid,ts,\
			tbgyn,tauyn,v0yn,widyn,tsyn,contyn) # co the sai thu tu 

	tb_tot = cfit.get_tb(xd, \
			tbgfit, taufit, v0fit, widfit, tsfit, cont_em)

	amg_oh['ts']    = tsfit ## DAU DAY
	amg_oh['tserr'] = tserr ## DAU DAY

	return amg_oh

## Gaussian fit - Step 0 ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def gfit_abs_line(cfit,x,bg,tau,v0,wid,cont,td,vrange,src,ell,bee,contoff):

	tfit, error, \
	tbgfit, taufit, v0fit, widfit, \
	tauerr, v0err, widerr, \
	cont,conterr,\
	cov, nloop, \
	nloopmax = cfit.abfit(x,td,vrange,cont,bg,tau,v0,wid)

	ngauss  = len(taufit)
	amg_oh  = {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
  		'tbg':0, 'tau':0.0, 'v0':0.0, 'wid':0.0, 'ts':0.0, \
  		'tbgerr':0, 'tauerr':0.0, 'v0err':0.0, 'widerr':0.0, 'tserr':0.0, \
  		'cont':0.0, 'conterr':0.0, 'contoff':0.0}
	
	dassign(amg_oh,'sname', src, ngauss)
	dassign(amg_oh,'ell', ell, ngauss)
	dassign(amg_oh,'bee', bee, ngauss)
	dassign(amg_oh,'gaussnr', list(range(ngauss)))
	dassign(amg_oh,'nrgauss', ngauss, ngauss)

	dassign(amg_oh,'tbg', tbgfit, ngauss)
	dassign(amg_oh,'tau', taufit)
	dassign(amg_oh,'v0', v0fit)
	dassign(amg_oh,'wid', widfit)

	dassign(amg_oh,'tauerr', tauerr.tolist())
	dassign(amg_oh,'v0err', v0err.tolist())
	dassign(amg_oh,'widerr', widerr)

	dassign(amg_oh,'cont', cont, ngauss)
	dassign(amg_oh,'conterr', conterr, ngauss)
	dassign(amg_oh,'contoff', contoff, ngauss)

	return amg_oh

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
		if(src != '3C105'):
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
		vrangid  = [xmin_id, xmax_id]
		vrangid  = [733,1023]

		## Background ##
		tbg1665    = 2.8+get_tb_408(ell[n], bee[n], inf408.tb_408)*(408./1665.)**2.8 # Tbg from 408MHz
		t_on1665   = td_ab
		t_off1665  = td_em
		tc1665     = cal_bg(xd,t_on1665)
		bg_off1665 = cal_bg(xd,t_off1665)
		trx1665    = bg_off1665 - tbg1665

		#Fitting, guess parameters and Fit multiple Gaussians for e(-tau)#
		guess = []
		for c in peak[src]:
			guess += peak[src][c]

		print n, '   ', src, ell[n], bee[n]

		cont       = 61.7   # Tc
		contoff    = 3.413776851462  # quan trong, thay doi la thay doi Tex rat nhieu
		tbaseline  = 52.215  # baseline of 1665 Emission line

		bg    = 0. # just 0. not to fit now
		tau   = [1.,1.]
		v0    = [8.2,10.2]
		wid   = [1.0,0.6]

		# amg_oh1 = gfit_abs_line(cfit,xd,bg,tau,v0,wid,cont,td_ab,vrangid,src,ell,bee,contoff)

		ta_fit, error, \
		xdat, ta_dat,\
		tbgfit, taufit, v0fit, widfit, \
		tauerr, v0err, widerr, \
		tcont,conterr,\
		cov, nloop, \
		nloopmax = cfit.abfit(xd,td_ab,vrangid,cont,bg,tau,v0,wid)

		ngauss  = len(taufit)
		amg_oh  = {'sname':'', 'ell':0.0, 'bee':0.0, 'gaussnr':0, 'nrgauss':1, \
	  		'tbg':0, 'tau':0.0, 'v0':0.0, 'wid':0.0, 'ts':0.0, \
	  		'tbgerr':0, 'tauerr':0.0, 'v0err':0.0, 'widerr':0.0, 'tserr':0.0, \
	  		'cont':0.0, 'conterr':0.0, 'contoff':0.0}
		
		dassign(amg_oh,'sname', src, ngauss)
		dassign(amg_oh,'ell', ell, ngauss)
		dassign(amg_oh,'bee', bee, ngauss)
		dassign(amg_oh,'gaussnr', list(range(ngauss)))
		dassign(amg_oh,'nrgauss', ngauss, ngauss)

		dassign(amg_oh,'tbg', tbgfit, ngauss)
		dassign(amg_oh,'tau', taufit)
		dassign(amg_oh,'v0', v0fit)
		dassign(amg_oh,'wid', widfit)

		dassign(amg_oh,'tauerr', tauerr.tolist())
		dassign(amg_oh,'v0err', v0err.tolist())
		dassign(amg_oh,'widerr', widerr)

		dassign(amg_oh,'cont', cont, ngauss)
		dassign(amg_oh,'conterr', conterr, ngauss)
		dassign(amg_oh,'contoff', contoff, ngauss)
		
		# fit with Tex
		ts      = [4.,7.]

		print tcont

		cont_em = contoff ## tbg=50./2

		tbg    = amg_oh['tbg'][0]
		tau    = amg_oh['tau']
		# tau    = [0.03, 0.05]
		v0     = amg_oh['v0']
		wid    = amg_oh['wid']
		npeaks = len(amg_oh['v0'])

		tdat   = td_em - tbaseline + contoff
		# tdat   = td_em

		te_fit, error, \
		xdat, te_dat,\
		toffset, tau, v0, wid, \
		cont,\
		cov, nloop, \
		nloopmax = cfit.emfit(xd,tdat,vrangid,tbaseline,tbg,tau,v0,wid,ts) # co the sai thu tu 

		etau_fit = ta_fit/tcont
		etau_dat = ta_dat/tcont

		ts_fit = (te_fit-contoff*etau_fit)/(1.-etau_fit)
		ts_dat = (te_dat-contoff*etau_dat)/(1.-etau_dat)

		plt.plot(xdat, ts_fit)
		plt.plot(xdat, ts_dat)
		plt.xlim(5.,12.)
		plt.ylim(-1.,30.)
		plt.show()

		# te_fit, error, \
		# toffset, tau, v0, wid, tsfit, \
		# tserr, \
		# cont,\
		# cov, nloop, \
		# nloopmax
		#cfit.emfit(xd,tdat,vrangid,cont_em,tbg,tau,v0,wid,ts) # co the sai thu tu 


		# print cont_em
		# amg_oh1['ts']    = tsfit ## DAU DAY
		# amg_oh1['tserr'] = tserr ## DAU DAY

		# print '1665 spin temps:'
		# print amg_oh1['ts']
		# print amg_oh1['tserr']

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
inf408 = readsav('../data/tb_408.sav') # l_cntr, b_cntr, tb_408
main(data, inf408)

sys.exit()