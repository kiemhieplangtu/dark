from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys, os

from scipy.io.idl import readsav


data = readsav('../data/makelines.sav')

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
			s   = s+stock[i]
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
        y = y - amp * np.exp( -((x - ctr)/wid)**2)
    return 1.-y


#==============
src_list = list(data.la.srcname)

for n in range(0,101) :

	sname = src_list[n]
	
	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	hi_f0 = data.la.cfr_bd0
	oh_f1 = data.la.cfr_bd1
	oh_f2 = data.la.cfr_bd2

	vlsr0 = data.la.vlsr_bd0
	vlsr1 = data.la.vlsr_bd1
	vlsr2 = data.la.vlsr_bd2

	em_avg0 = data.la.i_em_avg_bd0
	em_med0 = data.la.i_em_med_bd0
	ab_avg0 = data.la.i_abs_avg_bd0
	ab_med0 = data.la.i_abs_med_bd0

	em_avg1 = data.la.i_em_avg_bd1
	em_med1 = data.la.i_em_med_bd1
	ab_avg1 = data.la.i_abs_avg_bd1
	ab_med1 = data.la.i_abs_med_bd1

	em_avg2 = data.la.i_em_avg_bd2
	em_med2 = data.la.i_em_med_bd2
	ab_avg2 = data.la.i_abs_avg_bd2
	ab_med2 = data.la.i_abs_med_bd2

	em_avg2 = correct_ctrl_chnl(em_avg2)
	em_med2 = correct_ctrl_chnl(em_med2)
	ab_avg2 = correct_ctrl_chnl(ab_avg2)
	ab_med2 = correct_ctrl_chnl(ab_med2)

	# VLSR #
	x = vlsr2[n]

	# On-/Off source Tb and Background Continuum #
	t_on   = 0.5*ab_avg2[n]
	t_off  = 0.5*em_avg2[n]
	bg_on  = cal_bg(x,t_on)
	bg_off = cal_bg(x,t_off)

	# e(-tau) and tau #
	etau   = t_on/bg_on
	tau    = -np.log(etau)


	# plt.plot(velo, 2000.*tau)
	plt.plot(x, etau , 'r-')
	plt.title(sname+' '+str(oh_f2[n]))
	# plt.xlim(xmin,xmax)
	# plt.ylim(0.,200.)
	plt.show()