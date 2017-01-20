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

## Interpolate the VLSR=0 bin ##
 #
 # params list tb Temperature

 # return list tb Temperature
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def correct_ctrl_chnl(tb):
	for i in range(0,101):
		tb[i][1023] = (tb[i][1021] + tb[i][1022] + tb[i][1024] + tb[i][1025])/4.

	return tb

## Interpolate the VLSR=0 bin ##
 #
 # params list y y-data
 # params list x x-data
 # params float xmin X-min
 # params float xmax X-max

 # return float ymin & ymax
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def get_ylimits(y,x,xmin,xmax):
	max_id = get_vel_index(x, xmin)
	min_id = get_vel_index(x, xmax)
	ymin   = min(y[min_id:max_id])
	ymax   = max(y[min_id:max_id])

	return ymin, ymax

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
	src_info = vel.read()

	return src_info['src']

## Plot HI line and OH lines in same frame ##
 #
 # params dict data Data
 #
 # return Void
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def plot_hi_oh(data):
	src_list = list(data.la.srcname)

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

	em_avg0 = correct_ctrl_chnl(data.la.i_em_avg_bd0)
	em_med0 = correct_ctrl_chnl(data.la.i_em_med_bd0)
	ab_avg0 = correct_ctrl_chnl(data.la.i_abs_avg_bd0)
	ab_med0 = correct_ctrl_chnl(data.la.i_abs_med_bd0)

	em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	em_avg2 = correct_ctrl_chnl(correct_ctrl_chnl(em_avg2))
	em_med2 = correct_ctrl_chnl(correct_ctrl_chnl(em_med2))
	ab_avg2 = correct_ctrl_chnl(correct_ctrl_chnl(ab_avg2))
	ab_med2 = correct_ctrl_chnl(correct_ctrl_chnl(ab_med2))

	# Vrange infor #
	fname    = '../data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()

	# 26 src with no CO #
	fname    = '../data/26src_no_co.txt'
	cols     = ['idx','src']
	fmt      = ['i','s']

	src_no_co = restore(fname, 2, cols, fmt)
	s26info   = src_no_co.read()
	s26src    = s26info['src']

	for src in src_list:
		if (src != '3C123'):
			continue
		n = src_list.index(src)

		# VLSR #
		x0 = vlsr0[n]
		x1 = vlsr1[n]
		x2 = vlsr2[n]

		# On-/Off-source Tb and Background Continuum #
		t_on0   = 0.5*ab_avg0[n]
		t_off0  = 0.5*em_avg0[n]
		bg_on0  = cal_bg(x0,t_on0)
		bg_off0 = cal_bg(x0,t_off0)

		t_on1   = 0.5*ab_avg1[n]
		t_off1  = 0.5*em_avg1[n]
		bg_on1  = cal_bg(x1,t_on1)
		bg_off1 = cal_bg(x1,t_off1)

		t_on2   = 0.5*ab_avg2[n]
		t_off2  = 0.5*em_avg2[n]
		bg_on2  = cal_bg(x2,t_on2)
		bg_off2 = cal_bg(x2,t_off2)


		# e(-tau) and tau #
		etau0   = t_on0/bg_on0
		tau0    = -np.log(etau0)

		etau1   = t_on1/bg_on1
		tau1    = -np.log(etau1)

		etau2   = t_on2/bg_on2
		tau2    = -np.log(etau2)

		# Plot #
		xmin = -20.
		xmax = 30.

		ymin, ymax = get_ylimits(t_off1,x1,xmin,xmax)
		plt.subplot(4, 1, 1)
		plt.plot(x1, t_off1, 'k-')
		plt.title(src)
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		plt.ylabel('1665 Em (K)')
		plt.grid()

		plt.subplot(4, 1, 2)
		ymin, ymax = get_ylimits(etau1,x1,xmin,xmax)
		plt.plot(x1, etau1, 'k-')
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		plt.ylabel('1665 Abs (K)')
		plt.grid()

		plt.subplot(4, 1, 3)
		ymin, ymax = get_ylimits(t_off2,x2,xmin,xmax)
		plt.plot(x2, t_off2, 'b-')
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		plt.ylabel('1667 Em (K)')
		plt.title(src)
		plt.grid()

		plt.subplot(4, 1, 4)
		ymin, ymax = get_ylimits(etau2,x2,xmin,xmax)
		plt.plot(x2, etau2, 'b-')
		plt.xlabel('VLSR (km/s)')
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		plt.ylabel('1667 Abs (K)')
		plt.grid()

		plt.show()

#============== MAIN ==============#
data = readsav('../data/makelines.sav')
# hi_lines(data) # Plot HI line ONLY
plot_hi_oh(data) #Plot HI & OH lines

sys.exit()