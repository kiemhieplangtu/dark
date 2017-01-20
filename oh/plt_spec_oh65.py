import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt

from scipy.optimize      import curve_fit
from scipy.io.idl        import readsav
from numpy               import array
from restore             import restore
from plotting            import cplot

## Find the value of baseline ##
 #
 # params list vlsr VLSR
 # params list stock Stock-parameter

 # return float Value of baseline
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
 # params string src Source to be considered
 # params list srcs List of Sources
 # params dict data OH Data
 #
 # return Void
 #
 # version 11/2016 
 # author Nguyen Van Hiep ##
def get_spec_data(src, srcs, data):
	n     = srcs.index(src)
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
	x0 = vlsr0[n] ## HI
	x1 = vlsr1[n] ## OH 1665 MHz
	x2 = vlsr2[n] ## OH 1667 MHz

	# On-/Off-source Tb and Background Baseline #
	t_on0   = ab_avg0[n]        ## In fact: ab_avg = [Ton - Toff] = Tc*e(-tau)
	t_off0  = em_avg0[n]
	bg_on0  = cal_bg(x0,t_on0)  ## bg_on = Tc
	bg_off0 = cal_bg(x0,t_off0)

	t_on1   = ab_avg1[n]        ## In fact: ab_avg = [Ton - Toff] = Tc*e(-tau)
	t_off1  = em_avg1[n]
	bg_on1  = cal_bg(x1,t_on1)  ## bg_on = Tc
	bg_off1 = cal_bg(x1,t_off1)

	t_on2   = ab_avg2[n]        ## In fact: ab_avg = [Ton - Toff] = Tc*e(-tau)
	t_off2  = em_avg2[n]
	bg_on2  = cal_bg(x2,t_on2)  ## bg_on = Tc
	bg_off2 = cal_bg(x2,t_off2)

	return [x1,t_on1,t_off1, bg_on1, bg_off1]

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
	cols   = ['idx','src','avmin1','avmax1','avmin2','avmax2','evmin1','evmax1','evmin2','evmax2']
	fmt    = ['i',  's',  'f',      'f',    'f',      'f',     'f',    'f',     'f',      'f']
	dat    = restore(fname, 2, cols, fmt)
	vel    = dat.read(asarray=False)

	avmin1 = vel_info['avmin1']
	avmax1 = vel_info['avmax1']
	avmin2 = vel_info['avmin2']
	avmax2 = vel_info['avmax2']

	evmin1 = vel_info['evmin1']
	evmax1 = vel_info['evmax1']
	evmin2 = vel_info['evmin2']
	evmax2 = vel_info['evmax2']

	return avmin1[n],avmax1[n],avmin2[n],avmax2[n],evmin1[n],evmax1[n],evmin2[n],evmax2[n]

## Plot OH 1665 spectra - Emission and Absorption ##
 #
 # params dict data OH Data
 #
 # return Void
 #
 # version 11/2016 
 # author Nguyen Van Hiep ##
def plot_1665(data):
	srcs = list(data.srcname)
	print len(srcs), 'Sources:'
	print srcs

	src = raw_input('Choose a Source: ')
	if (src not in srcs):
		print 'Source not in the list. Please try again!'
		sys.exit()

	## Get infor VLSR, Ton, Toff, Tc, Bg_off
	v,ton,toff,bgon,bgoff = get_spec_data(src, srcs, data)

	# e(-tau) and tau #
	etau1   = ton/bgon
	tau1    = -np.log(etau1)

	# Plot #
	fig    = cplot()
	trace1 = fig.lines(v,etau1,label='OH 1665MHz',
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
	              xaxis  = dict(label='$V_{lsr} (km/s)$',tick_size=18,fontsize=35,xlim=[-30.,30.]),
	              yaxis  = dict(label=r'$e^{-\tau}$',tick_size=18,fontsize=35),
	              text   = [dict(loc=[-0.5,0.4],text='',color='blue',fontsize=17),
	              			dict(loc=[-0.5,0.31],text='',color='red',fontsize=19)
	              		   ],
	             )

	fig.iplot(data,layout,fullscreen=False)

	## Plot Emission line ##
	trace1 = fig.lines(v,toff,label='OH 1665MHz',
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
	              xaxis  = dict(label='$V_{lsr} (km/s)$',tick_size=18,fontsize=35,xlim=[-30.,30.]),
	              yaxis  = dict(label=r'$T(K)$',tick_size=18,fontsize=35),
	              text   = [dict(loc=[-0.5,0.4],text='',color='blue',fontsize=17),
	              			dict(loc=[-0.5,0.31],text='',color='red',fontsize=19)
	              		   ],
	             )

	fig.iplot(data,layout,fullscreen=False)

#============== MAIN ==============#
data = readsav('data/makelines.sav')
plot_1665(data.la)

sys.exit()