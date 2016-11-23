import sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

# Plot the correlation between Factors and log(N(HI)) only from Planck #
#
# params dict data N(HI) of 26 sources
#
# return Void
# 
# Author Van Hiep
##
def plot_planck_factor_vs_nhi(data):
	nhi_hl      = data['nhi_hl']
	nh_pl       = data['nh_pl']
	nhi_fk      = data['nhi_fk']

	err_pl      = data['err_pl']
	err_hl      = data['err_hl']
	err_fk      = data['err_fk']

	diff_fk_hl  = data['fk_hl']
	diff_pl_hl  = data['pl_hl']

	oh          = data['oh']
	src         = data['src']
	idx         = data['idx']

	# Error bar for x-axis and y-axis
	xerr_hl = np.array(err_hl)/np.array(nhi_hl)/np.log(10.0)
	yerr_pl = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nh_pl, err_pl)
	yerr_fk = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nhi_fk, err_fk)

	for k in range(0, len(nhi_hl)):
		nhi_hl[k] = np.log10(nhi_hl[k])

	# Fit and Plot #
	m,b = np.polyfit(nhi_hl,diff_pl_hl,1)
	n,c = np.polyfit(nhi_hl,diff_fk_hl,1)

	m,b,ea,eb = linear_fit(nhi_hl,diff_pl_hl)

	for i in range(len(nhi_hl)):
		print i, src[i], nhi_hl[i], xerr_hl[i], diff_pl_hl[i], yerr_pl[i], oh[i]

	fig    = cplot(filename='mmm')	
	trace2 = fig.error_bar(nhi_hl,diff_pl_hl,label='Factor $f = N_{H}$/$N_{HI}$',
			err_x=dict(val=xerr_hl),
			err_y=dict(val=yerr_pl),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)
	trace1 = fig.lines(nhi_hl,m*np.array(nhi_hl)+b,label='Best linear fit',
			prop=dict(color='b',
				      linewidth=3,
				      linestyle='solid',
				      marker='o',
				      markerfacecolor='b',
				      markersize=0
			))

	n = round(n,2)
	c = round(c,2)

	m = round(m,2)
	b = round(b,2)
	ea = round(ea,2)
	eb = round(eb,2)

	data   = [trace2,trace1]
	layout = dict(title  = 'Correlation between Factor $f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)',tick_size=18,fontsize=35,xlim=[0.,1.6]),
	              yaxis  = dict(label='$Factor f = N_{H}$/$N_{HI}$',tick_size=18,fontsize=35,ylim=[0.,3.]),
	              text   = [dict(loc=[0.2,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              			dict(loc=[0.2,0.3],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              		   ],
	             )

	for i,j in zip(idx,diff_pl_hl):
		if (oh[i] > 0) :
			fig.annotation(dict(label='('+str(src[i])+')',
				x=nhi_hl[i],y=j,
				fontsize=18,
				xycoords='data',
				xoff=-50.,yoff=30.,
				textcoords='offset points',
				arrowprops=dict(arrowstyle="->") )
				)

	fig.iplot(data,layout)

# Plot the correlation between Factors and log(N(HI)_CNM) only from Planck #
#
# params dict data N(HI) of 26 sources
#
# return Void
# 
# Author Van Hiep
##
def plot_planck_factor_vs_nhi_cnm(data):
	## Read N(HI)_fukui N(H)_Planck N(HI)_Heiles #
	cols     = ['idx','src','nhi_hl','wnm','cnm','er_hl','er_cnm']
	fmt      = ['i','s','f','f','f','f', 'f']
	dat      = restore('result/cnm_uncertainties_arcb.txt', 3, cols, fmt)
	cnm_inf  = dat.read() # the Uncertainties of CNM component

	nhi_cnm     = data['cnm']
	nhi_hl      = data['nhi_hl']
	nh_pl       = data['nh_pl']
	nhi_fk      = data['nhi_fk']

	err_pl      = data['err_pl']
	err_hl      = data['err_hl']
	err_fk      = data['err_fk']
	err_cnm     = cnm_inf['er_cnm']

	diff_fk_hl  = data['fk_hl']
	diff_pl_hl  = data['pl_hl']

	oh          = data['oh']
	src         = data['src']
	idx         = data['idx']

	dif_pl_cnm = list_div(nh_pl, nhi_cnm)

	# Error bar for x-axis and y-axis
	xerr_hl = np.array(err_cnm)/np.array(nhi_cnm)/np.log(10.0) # becaus x-axis = log10(N(HI))
	yerr_pl = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nh_pl, err_pl)

	for k in range(0, len(nhi_cnm)):
		nhi_cnm[k] = np.log10(nhi_cnm[k])

	# Fit and Plot #
	m,b,ea,eb = linear_fit(nhi_cnm,diff_pl_hl)

	# plt.errorbar(nhi_cnm, diff_pl_hl, yerr=yerr_pl, xerr=xerr_hl, fmt='ro', label='Factor $f = N_{H}$/$N_{HI}$')
	# plt.plot(nhi_cnm, m*np.array(nhi_cnm) + b, 'b-', label='Best linear fit', linewidth=2)

	fig    = cplot()
	trace2 = fig.error_bar(nhi_cnm,diff_pl_hl,label='Factor $f = N_{H}$/$N_{HI}$',
			err_x=dict(val=xerr_hl),
			err_y=dict(val=yerr_pl),			
			prop=dict(fmt='ro',
				markersize=8,
				markeredgecolor='b',
				markeredgewidth=1)
			)
	trace1 = fig.lines(nhi_cnm,m*np.array(nhi_cnm)+b,label='Best linear fit',
			prop=dict(color='b',
				      linewidth=3,
				      linestyle='solid',
				      marker='o',
				      markerfacecolor='b',
				      markersize=0
			))

	m = round(m,2)
	b = round(b,2)
	ea = round(ea,2)
	eb = round(eb,2)

	data   = [trace2,trace1]
	layout = dict(title  = 'Correlation between Factor $f = N_{H}$/$N_{HI}$ and CNM column density $N_{CNM}$ \nalong 26 lines-of-sight without the presence of CO line',
				  title_fontsize=30,
	              grid   = True,
	              legend = dict(loc='upper left', fontsize=18),
	              xaxis  = dict(label='log$_{10}$($N_{CNM}/10^{20}$ cm$^{-2}$)',tick_size=18,fontsize=35,xlim=[-0.6,1.28]),
	              yaxis  = dict(label='$Factor f = N_{H}$/$N_{HI}$',tick_size=18,fontsize=35,ylim=[0.,3.]),
	              text   = [dict(loc=[-0.5,0.4],text='a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb),color='blue',fontsize=17),
	              			dict(loc=[-0.5,0.31],text='(Available sources with the presence of OH line are shown)',color='red',fontsize=19)
	              		   ],
	             )

	for i,j in zip(idx,diff_pl_hl):
		if (oh[i] > 0) :
			fig.annotation(dict(label='('+str(src[i])+')',
				x=nhi_cnm[i],
				y=j,
				fontsize=12,
				xycoords='data',
				xoff=-50.,yoff=30.,
				textcoords='offset points',
				arrowprops=dict(arrowstyle="->") )
				)

	fig.iplot(data,layout)

# Calculate the uncertainties of factors #
#
# params 
# params 
# params 
#
# return list of factor_uncertainties
# 
# Author Van Hiep
##
def uncertainty_of_factors(factr, nhi_hl, err_hl, nhi, err):
	d1 = np.array(err_hl)/np.array(nhi_hl)
	d1 = np.square(d1)

	d2 = np.array(err)/np.array(nhi)
	d2 = np.square(d2)

	d  = np.sqrt(d1+d2)*np.array(factr)

	return d.tolist()

# plot N(HI)/N(HI)_thin vs N(HI)_thin #
#
# params dict data data to plot
#
# return void
# 
# Author Van Hiep
##
def plot_factor_vs_nhi(data):

	fact       = []
	lognhi     = []
	nhi        = data['nhi_thin'] # Optically-thin assumption
	nhi_heiles = data['nhi_hl'] # N(HI) from Carl Heiles Paper

	for i in range(0, len(nhi)):
		temp   = round(nhi_heiles[i]/nhi[i], 3)
		lognhi.append(np.log10(nhi[i]))
		fact.append(temp)

	# Fit and Plot #
	a, b = np.polyfit(lognhi,fact,1)

	plt.plot(lognhi, fact, 'r.', label='Factor = N(HI)/N(HI)_thin')
	plt.plot(lognhi, a*np.array(lognhi) + b, 'k-', label='')

	a = round(a, 2)
	b = round(b, 2)

	plt.xlabel('Factor')
	plt.ylabel('log(N(HI)_opt.thin).1e20')
	plt.title('Correlation between N(HI)_Heiles and Optically thin N(HI)')
	plt.grid(True)

	plt.text(0.21, 1.31, 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=12)
	plt.text(0.21, 0.92, 'a = '+str(a)+'  b = '+str(b), color='blue', fontsize=12)

	plt.legend(loc='upper right')
	plt.show()

# linear fit #
#
# params x list x-data
# params y list y-data
#
# return fit parameters and errors
# 
# Author Van Hiep
##
def linear_fit(x,y):
	sxy = 0.
	sx  = 0.
	sy  = 0.
	sx2 = 0.
	n   = len(x)
	for i in range(0,n) :
		sxy = sxy + x[i]*y[i]
		sx  = sx + x[i]
		sy  = sy + y[i]
		sx2 = sx2 + x[i]**2

	denom = (n*sx2 - sx**2)
	a = (n*sxy - sx*sy)/denom
	b = (sx2*sy - sx*sxy)/denom

	t    = n*sx2 - sx**2
	er_a = np.sqrt(n/t) 
	er_b = np.sqrt(sx2/t) 

	chi2 = 0.
	for i in range(0,n) :
		chi2 = chi2 + (y[i]-a*x[i]-b)**2

	chi2 = np.sqrt(chi2/(n-2))
	er_a = chi2*er_a
	er_b = chi2*er_b

	return a,b,er_a,er_b

## Divide two lists #
 #
 # params x list 
 # params y list 
 #
 # return x/y
 # 
 # Author Van Hiep
 ##
def list_div(x,y):
	ret = []
	for i in range(0, len(x)):
		ret.append(x[i]/y[i])

	return ret

#================= MAIN ========================#

## Read N(HI)_fukui N(H)_Planck N(HI)_Heiles #
 # Author Van Hiep
 ##
cols = ['idx','src','nhi_fk','err_fk','nh_pl','err_pl','nhi_hl','er_hl','err_hl','wnm','cnm','fk_hl','pl_hl','oh','nhi_thin']
fmt  = ['i','s','f','f','f','f','f','f','f','f','f','f','f','f','f']
data = restore('result/26src_no_co_nhi_and_uncertainties_full.txt', 4, cols, fmt)
dat  = data.read()

plot_planck_factor_vs_nhi(dat)
# plot_planck_factor_vs_nhi_cnm(dat)