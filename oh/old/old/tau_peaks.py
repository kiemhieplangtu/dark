import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot

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
		tb[i][1023] = (tb[i][1021] + tb[i][1022] + tb[i][1024] + tb[i][1025])/4.

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
 # params float bg_off  Off-source Background_continuum
 # params list popt     Fit-result parameters 
 #
 # return List tex      Excitation temperature
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tex_from_etau(etau_fit,velo,tb_off,bg_off,popt):
	tex  = []

	# Find the vel-range of the peak #
	peak = cal_vel_range_fwhm(popt)

	# Compute the Excitation Temperature #
	s     = 0.
	count = 0
	for i in range(0, len(etau_fit)):
		if (etau_fit[i] != 1.):
			ts = (tb_off[i]-bg_off*etau_fit[i])/(1.-etau_fit[i])
			tex.append(ts)
		else:
			tex.append(0.)

		for k in range(0,len(peak),2):
			vmin = peak[0+k]
			vmax = peak[1+k]

			if ((velo[i]>vmin) and (velo[i]<vmax)) :
				s     = s + tex[i]
				count = count + 1

	s = s/count
	return tex,s

## Create multiple Gaussian functions ##
 #
 # params string fname Filename #
 # return list List of MM_Survey Sources
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_mm_survey_src(fname = ''):
	cols     = ['src']
	fmt      = ['s']
	inf      = restore(fname, 0, cols, fmt)
	src_info = inf.read()

	return src_info['src']

## Compute Tau peaks for 1667 line ##
 #
 # params dict data Data
 #
 # return list Tau peaks for 1667
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tau_peak_1667(data):
 	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	oh_f2 = data.la.cfr_bd2
	vlsr2 = data.la.vlsr_bd2

	em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	# Gaussian peaks' info - estimate values of Amp, tau0 and Width #
	peak = peak_info('data/gauss_1667_peaks.txt')

	# Vrange infor #
	fname    = 'data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']

	# src of MM Survey #
	mm_src   = get_mm_survey_src('data/mm_survey_src.txt')

	tau_peak = []
	for src in src_list:
		n    = src_list.index(src)
		xmin = vmin[n]
		xmax = vmax[n]

		if ((xmin == 0. and xmax == 0.) or (src not in mm_src)):
			continue

		# VLSR #
		x = vlsr2[n]

		# On-/Off-source Tb and Background Continuum #
		t_on   = 0.5*ab_avg2[n]
		t_off  = 0.5*em_avg2[n]
		bg_on  = cal_bg(x,t_on)
		bg_off = cal_bg(x,t_off)

		# e(-tau) and tau #
		etau   = t_on/bg_on
		tau    = -np.log(etau)

		# Get Index of the Velocity-range #
		xmax_id  = get_vel_index(x, xmin)   # xmax_index
		xmin_id  = get_vel_index(x, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id          # Total number of bins

		# Get other quantities within the Vel-range
		velo    = []
		e_tau   = []
		tb_off  = []
		for i in range(0,num_chnl):
			velo.append(x[xmin_id+i])
			e_tau.append(etau[xmin_id+i])
			tb_off.append(t_off[xmin_id+i])

		#Fitting, guess parameters and Fit multiple Gaussians for e(-tau)#
		guess = []
		for c in peak[src]:
			guess += peak[src][c]

		popt, pcov = curve_fit(func, velo, e_tau, p0=guess) # Do the fit
		# print popt

		vpeak = []
		for j in range(0, len(popt), 3):
			ctrl_vel = popt[1+j]
			vpeak    += [ctrl_vel]    # vertical lines at 2*sigma

		etau_fit = func(velo, *popt) # e(-tau) of the fit
		etau_pk  = func(vpeak, *popt) # e(-tau) of the fit
		tau      = -np.log(etau_fit) # tau

		for j in range(0, len(popt), 3):
			ctrl_vel = popt[1+j]
			etau_pk  = func(ctrl_vel, *popt)
			tau_pk   = -np.log(etau_pk)
			tau_peak += [tau_pk]

	return tau_peak	
	
## Compute Tau peaks for 1665 line ##
 #
 # params dict data Data
 #
 # return list Tau peaks for 1665
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tau_peak_1665(data):
 	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	oh_f  = data.la.cfr_bd1
	vlsr  = data.la.vlsr_bd1

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

	# src of MM Survey #
	mm_src   = get_mm_survey_src('data/mm_survey_src.txt')

	tau_peak = []
	for src in src_list:
		n    = src_list.index(src)
		xmin = vmin[n]
		xmax = vmax[n]

		if ((xmin == 0. and xmax == 0.) or (src not in mm_src)):
			continue

		# VLSR #
		x = vlsr[n]

		# On-/Off-source Tb and Background Continuum #
		t_on   = 0.5*ab_avg[n]
		t_off  = 0.5*em_avg[n]
		bg_on  = cal_bg(x,t_on)
		bg_off = cal_bg(x,t_off)

		# e(-tau) and tau #
		etau   = t_on/bg_on
		tau    = -np.log(etau)

		# Get Index of the Velocity-range #
		xmax_id  = get_vel_index(x, xmin)   # xmax_index
		xmin_id  = get_vel_index(x, xmax)   # xmin_index
		num_chnl = xmax_id-xmin_id          # Total number of bins

		# Get other quantities within the Vel-range
		velo    = []
		e_tau   = []
		tb_off  = []
		for i in range(0,num_chnl):
			velo.append(x[xmin_id+i])
			e_tau.append(etau[xmin_id+i])
			tb_off.append(t_off[xmin_id+i])

		#Fitting, guess parameters and Fit multiple Gaussians for e(-tau)#
		guess = []
		for c in peak[src]:
			guess += peak[src][c]

		popt, pcov = curve_fit(func, velo, e_tau, p0=guess) # Do the fit
		# print popt

		vpeak = []
		for j in range(0, len(popt), 3):
			ctrl_vel = popt[1+j]
			vpeak    += [ctrl_vel]    # vertical lines at 2*sigma

		etau_fit = func(velo, *popt) # e(-tau) of the fit
		etau_pk  = func(vpeak, *popt) # e(-tau) of the fit
		tau      = -np.log(etau_fit) # tau

		for j in range(0, len(popt), 3):
			ctrl_vel = popt[1+j]
			etau_pk  = func(ctrl_vel, *popt)
			tau_pk   = -np.log(etau_pk)
			tau_peak += [9.*tau_pk/5.]

		# 	if((tau_pk > 0.28) or (tau_pk < 0.)):
		# 		print src, 9.*tau_pk/5.

		# print '======'
		# print '======'
	return tau_peak

#============== MAIN ==============#
data = readsav('data/makelines.sav')
peak_tau_1667 = cal_tau_peak_1667(data)
peak_tau_1665 = cal_tau_peak_1665(data)

# Plot histogram #
size = 0.01
fig  = cplot()
trace1 = fig.hist(np.asarray(peak_tau_1667),label='1667 tau peaks',autobinx=False,
                  xbins=dict(start=0.0, end=0.3, size=size),
                  opacity=1.0,
                  marker=dict(
                    color = 'b',                    
                    )
                 )
trace2 = fig.hist(np.asarray(peak_tau_1665),label='1665 tau peaks',autobinx=False,
                  xbins=dict(start=0.0, end=0.3, size=size),
                  opacity=1.0,
                  marker=dict(
                    color = 'g',                    
                    )
                 )
data   = [trace2,trace1]

layout = dict(title  = 'Histogram of Tau peaks',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper right', fontsize=18),
	              xaxis  = dict(label='Tau peak values',tick_size=18,fontsize=35),
	              yaxis  = dict(label='Counts',tick_size=18,fontsize=35),
	              # text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              # 			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              # 		   ],
	             )
fig.iplot(data,layout)

sys.exit()