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
        y = y - amp * np.exp( -((x - ctr)/wid)**2)
    return 1.-y

## Gaussian fit for 1667 line ##
 #
 # params dict data Data
 #
 # return fit-results
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def hi_lines(data):
 	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	oh_f = data.la.cfr_bd0
	vlsr = data.la.vlsr_bd0

	em_avg = correct_ctrl_chnl(data.la.i_em_avg_bd0)
	em_med = correct_ctrl_chnl(data.la.i_em_med_bd0)
	ab_avg = correct_ctrl_chnl(data.la.i_abs_avg_bd0)
	ab_med = correct_ctrl_chnl(data.la.i_abs_med_bd0)

	# Peaks of considered sources #
	peak = {}
	peak['3C18']  = {'0':[-0.013,-8.04, 0.7]}
	peak['3C109'] = {'0':[-0.0037, 9.22, 0.52],'1':[-0.0056, 10.06, 0.63]}
	peak['3C409'] = {'0':[-0.0292, 15.35,0.62]}
	peak['3C132'] = {'0':[-0.0057,7.794,0.518]}
	peak['P0531+19'] = {'0':[-0.0011,1.88,0.9],'1':[-0.00095,5.047,1.0]}

	# Vrange infor #
	# fname    = '../data/vel_range.txt'
	# cols     = ['idx','src','vmin','vmax']
	# fmt      = ['i','s','f','f']
	# vel      = restore(fname, 2, cols, fmt)
	# vel_info = vel.read()
	# vmin    = vel_info['vmin']
	# vmax    = vel_info['vmax']

	# 26 src with no CO #
	# fname    = '../data/26src_no_co.txt'
	# cols     = ['idx','src']
	# fmt      = ['i','s']
	# src_no_co = restore(fname, 2, cols, fmt)
	# s26info   = src_no_co.read()
	# s26src    = s26info['src']

	for src in src_list:
		n = src_list.index(src)

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

		# Plot #
		fig    = cplot()
		trace1 = fig.lines(x,etau,label='exp(-tau)',
				prop=dict(color='r',
					      linewidth=1,
					      linestyle='solid',
					      marker='o',
					      markerfacecolor='b',
					      markersize=0
				))

		data   = [trace1]
		layout = dict(title  = 'HI 21 cm - '+str(n)+'th: '+src+'  exp(-tau)',
					  title_fontsize=30,
		              grid   = True,
		              legend = dict(loc='upper left', fontsize=18),
		              xaxis  = dict(label='vlsr (km/s)',tick_size=18,fontsize=35,xlim=[-60.,60.]),
		              yaxis  = dict(label='exp(-tau)',tick_size=18,fontsize=35),
		              text   = [dict(loc=[-0.5,0.4],text='',color='blue',fontsize=17),
		              			dict(loc=[-0.5,0.31],text='',color='red',fontsize=19)
		              		   ],
		             )

		fig.iplot(data,layout)


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

	# Peaks of considered sources #
	peak = {}
	peak['3C18']  = {'0':[-0.014,-7.8, 0.5]}
	peak['3C109'] = {'0':[-0.03, 9.2, 0.5],'1':[-0.055, 10.5, 0.6]}
	peak['3C409'] = {'0':[-0.03	,15.5,0.5]}
	peak['3C132'] = {'0':[-0.07,7.7,0.5]}
	peak['P0531+19'] = {'0':[-0.002,1.8,0.7],'1':[-0.002,5.3,0.8]}


	# Vrange infor #
	fname    = '../data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin    = vel_info['vmin']
	vmax    = vel_info['vmax']

	# 26 src with no CO #
	fname    = '../data/26src_no_co.txt'
	cols     = ['idx','src']
	fmt      = ['i','s']

	src_no_co = restore(fname, 2, cols, fmt)
	s26info   = src_no_co.read()
	s26src    = s26info['src']

	for src in src_list:
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

		etau0  = (etau0-1.)*0.025+0.06
		etau1  = (etau1-1.)+0.03
		etau2  = (etau2-1.)

		# Plot #
		fig    = cplot()
		trace0 = fig.lines(x0,etau0,label='exp(-tau) - HI',
				prop=dict(color='b',
					      linewidth=1,
					      linestyle='solid',
					      marker='o',
					      markerfacecolor='b',
					      markersize=0
				))

		trace1 = fig.lines(x1,etau1,label='exp(-tau) - OH(1665)',
				prop=dict(color='r',
					      linewidth=1,
					      linestyle='solid',
					      marker='o',
					      markerfacecolor='b',
					      markersize=0
				))

		trace2 = fig.lines(x2,etau2,label='exp(-tau) - OH(1667)',
				prop=dict(color='purple',
					      linewidth=1,
					      linestyle='solid',
					      marker='o',
					      markerfacecolor='b',
					      markersize=0
				))

		data   = [trace0, trace1, trace2]
		layout = dict(title  = 'All lines - '+str(n)+'th: '+src+'  exp(-tau)',
					  title_fontsize=30,
		              grid   = True,
		              legend = dict(loc='upper left', fontsize=18),
		              xaxis  = dict(label='vlsr (km/s)',tick_size=18,fontsize=35,xlim=[-60.,60.]),
		              yaxis  = dict(label='exp(-tau)',tick_size=18,fontsize=35),
		              text   = [dict(loc=[-0.5,0.4],text='',color='blue',fontsize=17),
		              			dict(loc=[-0.5,0.31],text='',color='red',fontsize=19)
		              		   ],
		             )

		fig.iplot(data,layout)

#============== MAIN ==============#
data = readsav('../data/makelines.sav')
# hi_lines(data)
plot_hi_oh(data)

sys.exit()






























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

# Peaks of considered sources #
peak = {}
peak['3C18']  = {'0':[-0.014,-7.8, 0.5]}
peak['3C109'] = {'0':[-0.03, 9.2, 0.5],'1':[-0.055, 10.5, 0.6]}
peak['3C409'] = {'0':[-0.03	,15.5,0.5]}
peak['3C132'] = {'0':[-0.07,7.7,0.5]}
peak['P0531+19'] = {'0':[-0.002,1.8,0.7],'1':[-0.002,5.3,0.8]}


# Vrange infor #
fname    = '../data/vel_range.txt'
cols     = ['idx','src','vmin','vmax']
fmt      = ['i','s','f','f']
vel      = restore(fname, 2, cols, fmt)
vel_info = vel.read()
vmin    = vel_info['vmin']
vmax    = vel_info['vmax']

# 26 src with no CO #
fname    = '../data/26src_no_co.txt'
cols     = ['idx','src']
fmt      = ['i','s']

src_no_co = restore(fname, 2, cols, fmt)
s26info   = src_no_co.read()
s26src    = s26info['src']

for src in s26src:
	if(src in src_list):
		n = src_list.index(src)
		xmin = vmin[n]
		xmax = vmax[n]

		if (xmin == 0. and xmax == 0.):
			continue
		print n,src,xmin,xmax

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
		xmax_id  = get_vel_index(x, xmin)
		xmin_id  = get_vel_index(x, xmax)
		num_chnl = xmax_id-xmin_id # Total number of bins

		# Get other quantities within the Vel-range
		velo    = []
		e_tau   = []
		tb_off  = []
		tex     = []
		for i in range(0,num_chnl):
			velo.append(x[xmin_id+i])
			e_tau.append(etau[xmin_id+i])
			tb_off.append(t_off[xmin_id+i])

		#Fitting, guess parameters and Fit multiple Gaussians for e(-tau)#
		guess = []
		for c in peak[src]:
			print c, peak[src][c]
			guess += peak[src][c]

		popt, pcov = curve_fit(func, velo, e_tau, p0=guess) # Do the fit
		print popt

		etau_fit = func(velo, *popt) # e(-tau) of the fit
		tau      = -np.log(etau_fit) # tau

		#Compute the Excitation Temperature #
		for i in range(0, len(etau_fit)):
			if (etau_fit[i] != 1.):
				ts = (tb_off[i]-bg_off*etau_fit[i])/(1.-etau_fit[i])
				# if(np.abs(ts) > 200.):
				# 	ts = 0.
				tex.append(ts)
			else:
				tex.append(0.)


		# plt.plot(velo, e_tau)
		# plt.plot(velo, etau_fit, 'r-')
		# plt.xlim(xmin,xmax)
		# plt.ylim(0.,200.)

		plt.plot(velo, 3000.*tau)
		plt.plot(velo, tex , 'r-')
		plt.xlim(xmin,xmax)
		plt.ylim(0.,200.)

		plt.title(str(n)+': '+src)
		plt.show()