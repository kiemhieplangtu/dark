import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import copy

from scipy.io.idl        import readsav
from numpy               import array
from restore             import restore
from plotting            import cplot
from mpfit               import mpfit

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
 	return slope, bsline, sigma

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

## Prepare data for OH1665 ##
 #
 # params int n Index of source
 # params dict data OH Data
 # params list vlims Velo-ranges to cal. baseline bckground
 #
 # return Void
 #
 # version 11/2016 
 # author Nguyen Van Hiep ##
def get_spec_data(n, data,vlims):
	ra50  = data.ra1950
	dec50 = data.dec1950
	ell   = data.ell
	bee   = data.bee

	hi_f0 = data.cfr_bd0
	oh_f1 = data.cfr_bd1
	oh_f2 = data.cfr_bd2

	vlsr0 = data.vlsr_bd0
	vlsr1 = data.vlsr_bd1
	vlsr2 = data.vlsr_bd2

	em_avg0 = correct_ctrl_chnl(data.i_em_avg_bd0)
	em_med0 = correct_ctrl_chnl(data.i_em_med_bd0)
	ab_avg0 = correct_ctrl_chnl(data.i_abs_avg_bd0)
	ab_med0 = correct_ctrl_chnl(data.i_abs_med_bd0)

	em_avg1 = correct_ctrl_chnl(data.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.i_abs_med_bd1)

	em_avg2 = correct_ctrl_chnl(data.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.i_abs_med_bd2)

	# VLSR #
	x0 = vlsr0[n] ## HI 1420 MHz
	x1 = vlsr1[n] ## OH 1665 MHz
	x2 = vlsr2[n] ## OH 1667 MHz

	# avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2 = read_vrange_to_cal_bg(n)

	# On-/Off-source Tb and Background Continuum #
	t_on0   = ab_avg0[n]        ## In fact: ab_avg = [Ton - Toff] = Tc*e(-tau)
	t_off0  = em_avg0[n]
	slope0,bg_on0,er0  = baseline_from_linear_fit(x0, t_on0, vlims[0],vlims[1],vlims[2],vlims[3],fit=False) ## bg_on = Tc
	slope0,bg_off0,er0 = baseline_from_linear_fit(x0, t_off0, vlims[4],vlims[5],vlims[6],vlims[7],fit=False)

	t_on1   = ab_avg1[n]        ## In fact: ab_avg = [Ton - Toff] = Tc*e(-tau)
	t_off1  = em_avg1[n]
	slope1,bg_on1,er1on  = baseline_from_linear_fit(x1, t_on1, vlims[0],vlims[1],vlims[2],vlims[3],fit=False) ## bg_on = Tc
	slope1,bg_off1,er1off = baseline_from_linear_fit(x1, t_off1, vlims[4],vlims[5],vlims[6],vlims[7],fit=False)

	t_on2   = ab_avg2[n]        ## In fact: ab_avg = [Ton - Toff] = Tc*e(-tau)
	t_off2  = em_avg2[n]
	slope2,bg_on2,er2  = baseline_from_linear_fit(x2, t_on2, vlims[0],vlims[1],vlims[2],vlims[3],fit=False) ## bg_on = Tc
	slope2,bg_off2,er2 = baseline_from_linear_fit(x2, t_off2, vlims[4],vlims[5],vlims[6],vlims[7],fit=False)

	return [x1,t_on1,t_off1, bg_on1, bg_off1, er1on, er1off]

## Read vel-range to calculate background (Absorption & Emission lines) ##
 #
 # params int n Order of the Source
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_vrange_to_cal_bg(n,fname='../sub_data/bins_to_cal_bg.txt'):
	cols  = ['idx','src','avmin1','avmax1','avmin2','avmax2','evmin1','evmax1','evmin2','evmax2']
	fmt   = ['i',  's',  'f',      'f',    'f',      'f',     'f',    'f',     'f',      'f']
	dat   = restore(fname, 2, cols, fmt)
	vel   = dat.read(asarray=False)

	amin1 = vel['avmin1']
	amax1 = vel['avmax1']
	amin2 = vel['avmin2']
	amax2 = vel['avmax2']

	emin1 = vel['evmin1']
	emax1 = vel['evmax1']
	emin2 = vel['evmin2']
	emax2 = vel['evmax2']

	return amin1[n],amax1[n],amin2[n],amax2[n],emin1[n],emax1[n],emin2[n],emax2[n]

## Get Velo-range of the lines ##
 #
 # params int n Index/order of Source
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

	v  = np.asarray(v, dtype=np.float64)
	tb = np.asarray(tb, dtype=np.float64)

 	return np.std(tb)

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

## Read 1666MHz-continuum Tb from 21cm continuum Healpy map ##
 #
 # params string fname File-name
 # return dict ret1, ret2: 1666MHz-continuum Tb
 #
 # version 11/2016 
 # author Nguyen Van Hiep ##
def read_tbg_cont(fname='../../hi/stockert_fits/result/tbg_cont_from_1420.txt'):
	cols = ['idx','src','l','b', 'tb21', 'tbc1', 'tbc2', 'tb408']
	fmt  = ['i',   's', 'f','f',   'f',   'f',    'f',    'f'   ]
	src  = restore(fname, 5, cols, fmt)
	info = src.read()

	src  = info['src']
	gl   = info['l']
	gb   = info['b']
	tc1  = info['tbc1']
	tc2  = info['tbc2']

	ret1 = {}
	ret2 = {}
	for i in range(len(src)):
		ret1[src[i]] = tc1[i]
		ret2[src[i]] = tc2[i]

	return ret1, ret2

## Plot OH 1665 spectra - Emission and Absorption ##
 #
 # params dict data OH Data
 #
 # return Void
 #
 # version 11/2016 
 # author Nguyen Van Hiep ##
def plot_1665(data, bd=1):
	## Continuum Background of 1665MHz (derived from 21cm Continuum) ##
	tc1,tc2 = read_tbg_cont('../../hi/stockert_fits/result/tbg_cont_from_1420.txt')

	## Absorption and Emission data ##
	cst = 3.99757843817
	frq = 1665.402
	pfl = '../data/gauss_1665_peaks.txt'
	if(bd == 2):
		cst = 2.21841824609
		frq = 1667.359
		pfl = '../data/gauss_1667_peaks.txt'

	srcs    = list(data.srcname)
	src     = '3C18'
	n       = srcs.index(src)

	## Get infor VLSR, Ton, Toff, Tc, Bg_off
	vlims  = read_vrange_to_cal_bg(n)
	avmin1 = vlims[0]
	avmax1 = vlims[1]
	avmin2 = vlims[2]
	avmax2 = vlims[3]
	evmin1 = vlims[4]
	evmax1 = vlims[5]
	evmin2 = vlims[6]
	evmax2 = vlims[7]

	v,ton,toff,bgon,bgoff,eron,eroff = get_spec_data(n,data,vlims)
 	xmin, xmax                       = vel_range(n)
 	if (xmin == 0. and xmax == 0.): sys.exit() #continue

	# if(src=='3C123'):
	# 	tbg1665 = 26.

	## BACKGROUNDS and THEIR UNCERTAINTIES ##
	tc       = bgon
	tc_er    = eron
	bgoff    = bgoff
	bgoff_er = eroff
	tcont21  = tc1[src]  ## 2.8+tb[i]*(1420./1665.402)**2.8
	trx      = bgoff - tcont21

	## 1-SIGMA STANDARD DEVIATION OF Tabspt and Temmission ##
	tab_sigma  = tc_er
	tem_sigma  = bgoff_er
	trx_sigma  = bgoff_er

	## COMPUTE EXP(-TAU), TAU & 1-SIGMA STANDARD SEVIATION OF TAU ##
	etaud      = ton/tc
	taud       = -np.log(etaud)
	tau_sigma  = get_tb_sigma(v, taud, avmin1, avmax1, avmin2, avmax2)
	etau_sigma = np.abs(etaud)*np.sqrt( (tab_sigma/ton)**2 + (tc_er/tc)**2 )
	etau_sigma = get_tb_sigma(v, etaud, avmin1, avmax1, avmin2, avmax2)

	# VELOCITY-RANGE & INDEXES #
	xmax_id  = get_vel_index(v, xmin)   # xmax_index
	xmin_id  = get_vel_index(v, xmax)   # xmin_index
	num_chnl = xmax_id-xmin_id           # Total number of bins
	vrange   = [xmin_id, xmax_id]
	dv       = (xmax-xmin)/num_chnl

	## (FOR FUN) FIT ABSORPTION LINE FOR TAU, V0 and WIDTH ##
	guesspar,base_range = peak_info(src,pfl)
	lguess              = [tc] + guesspar
	x,etaufit,etau,\
	abp,abper,npar,\
	parbase,pname,parinfo = ab_fit(src,v,ton,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2)

	## PLOT FIT ##
	# colors = ['m','g','b','y','c','r','purple','b']
	# plt.plot(x,etau)
	# plt.plot(x,etaufit,'r-')
	# for i in range(2,len(abp),3):
	# 	plt.axvline(abp[i]-abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color=colors[(i-3)/4])
	# 	plt.axvline(abp[i]+abp[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color=colors[(i-3)/4])
	# plt.title('Absorption fit - ' + src)
	# plt.grid()
	# plt.show()
	## END: PLOT FIT ##

	## Tau and Width to cal. N(OH) ##
	tau_fit = []
	v0_fit  = []
	wid_fit = []
	for i in range(1,len(abp),3):
		tau_fit.append(abp[i])
		v0_fit.append(abp[i+1])
		wid_fit.append(abp[i+2])

	## CALCULATE Tex, CHOOSE etaufit OR etau ##
	t_on     = toff + ton # On-source Spectrum
	t_on     = t_on[xmin_id:xmax_id]
	xde      = v[xmin_id:xmax_id]
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
	tex     = (ton-trx-(tcont21 + tc)*etaue)/(1.-etaue)

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

	sys.exit()

	# e(-tau) and tau #
	etau1   = ton/bgon
	tau1    = -np.log(etau1)

	
	# Plot #
	fig    = cplot()
	trace1 = fig.lines(v,etau1,label='OH 1665',
			prop=dict(color='b',
				      linewidth=1,
				      linestyle='solid',
				      marker='o',
				      markerfacecolor='b',
				      markersize=0
			))

	data   = [trace1]
	layout = dict(title  = 'Absorption line: ' + src,
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='vlsr (km/s)',tick_size=18,fontsize=35,xlim=[-60.,60.]),
	              yaxis  = dict(label=r'$e^{-\tau}$',tick_size=18,fontsize=35),
	              text   = [dict(loc=[-0.5,0.4],text='',color='blue',fontsize=17),
	              			dict(loc=[-0.5,0.31],text='',color='red',fontsize=19)
	              		   ],
	             )

	fig.iplot(data,layout,fullscreen=False)

	## Plot Emission line ##
	trace1 = fig.lines(v,toff,label='OH 1665',
			prop=dict(color='b',
				      linewidth=1,
				      linestyle='solid',
				      marker='o',
				      markerfacecolor='b',
				      markersize=0
			))

	data   = [trace1]
	layout = dict(title  = 'Emission line: ' + src,
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='vlsr (km/s)',tick_size=18,fontsize=35,xlim=[-60.,60.]),
	              yaxis  = dict(label=r'$T(K)$',tick_size=18,fontsize=35),
	              text   = [dict(loc=[-0.5,0.4],text='',color='blue',fontsize=17),
	              			dict(loc=[-0.5,0.31],text='',color='red',fontsize=19)
	              		   ],
	             )

	fig.iplot(data,layout,fullscreen=False)

#============== MAIN ==============#
data = readsav('../data/makelines.sav')
plot_1665(data.la)

sys.exit()