import os, sys, shutil
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import pymc3             as pm
import copy

from numpy   import array
from restore import restore
from mpfit   import mpfit
import operator

## Linear ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
# def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
def myfunc(p, fjac=None, x=None, y=None, err=None):
	model  = p[0] * x + p[1]
	status = 0
	return [status, (y - model) / err]

## Read N(HI)_Heiles for 26 src without CO #
 # Then calculate the uncertainties for each component
 #
 # params string fname Filename
 #
 # return dict Gaussian-fit parameters of each component
 # including the uncertainties of each component
 # 
 # version 10/2016
 # Author Van Hiep ##
def read_nhi_fukui_nh_planck(fname = 'result/26src_no_co_nhi_and_uncertainties_full.txt'):
	cols  = ['idx','src','nhi_fk','err_fk','nh_pl','err_pl','nhi_hl','er_hl', 'err_hl', 'wnm', 'cnm', 'fk_hl', 'pl_hl', 'oh', 'nhi_thin']
	fmt   = ['i',   's',  'f',    'f',      'f',    'f',     'f',    'f',       'f',     'f',   'f',   'f',    'f',       'f',    'f']
	dat   = restore(fname, 4, cols, fmt)
	inf   = dat.read()
	return inf

## Read 23 src with low N(HI) #
 #
 # params string fname Filename
 #
 # return dict Infor
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_23src_lownhi(fname = 'result/23src_lownhi_nhi_and_uncertainties.txt'):
	cols  = ['idx','src','nhi_fk','err_fk','nh_pl','err_pl','nhi_hl','er_hl', 'err_hl', 'wnm', 'cnm', 'fk_hl', 'pl_hl', 'oh', 'nhi_thin']
	fmt   = ['i',   's',  'f',    'f',      'f',    'f',     'f',    'f',       'f',     'f',   'f',   'f',    'f',       'f',    'f']
	dat   = restore(fname, 4, cols, fmt)
	inf   = dat.read()
	return inf

## Read Uncertainties of CNM (& WNM) #
 #
 # params string fname Filename
 # return dict info of N(H) and N(HI)
 # 
 # version 10/2016
 # Author Van Hiep ##
def read_cnm_err(fname = '../hi/result/26src_no_co_cnm_uncertainties_arcb.txt'):
	cols  = ['idx','src','nhi','nhi_er','cnm','cnm_er','wnm','wnm_er']
	fmt   = ['i',   's',  'f', 'f',     'f',    'f',     'f',     'f']
	dat   = restore(fname, 2, cols, fmt)
	inf   = dat.read()
	return inf

## Plot the correlation between Factors and log(N(HI)) only from Planck #
 #
 # params dict data N(HI) of 26 sources without CO
 # params dict data N(HI) of 23 sources with Low NHI
 #
 # return Void
 # 
 # Author Van Hiep ##
def plot_planck_factor_vs_nhi(data, lownhi):
	## 26 wrc without CO
	nhi         = data['nhi_hl']
	nhi_hl      = nhi
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

	hi          = lownhi['nhi_hl']
	nh          = lownhi['nh_pl']
	err_hi      = lownhi['err_hl']
	err_nh      = lownhi['err_pl']

	print nhi
	# plt.plot(nhi,nh_pl, 'rd', label='data', ms=10)
	plt.errorbar(nhi,nh_pl,xerr=err_hl, yerr=err_pl, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data no CO')
	plt.errorbar(hi,nh,xerr=err_hi, yerr=err_nh, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data low $N_{HI}, N_{HI}<3.0\cdot10^{20} cm^{-2}$')
	plt.plot([0,30],[0,30], 'k--', label='$N_{H} = N_{HI}$')
	plt.title('Correlation between $N_{H}$ and $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line & 23 LOS with low $N_{HI}$', fontsize=30)
	plt.ylabel('$N_{H}[10^{20}$ cm$^{-2}]$', fontsize=35)
	plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$', fontsize=35)
	# plt.xlim(0, 1.6)
	# plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(15., 2., '(Available sources with the presence of OH are shown)', color='k', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(nhi[i], nh_pl[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=18,
	               )

	plt.show()

## Plot the correlation between Factors and log(N(HI)_CNM) only from Planck #
 #
 # params dict data N(HI) of 26 sources
 # return Void
 # 
 # Author Van Hiep ##
def plot_planck_factor_vs_nhi_cnm(data):

	cnm_inf     = read_cnm_err('../hi/result/26src_no_co_cnm_uncertainties_arcb.txt') # the Uncertainties of CNM component

	nhi_cnm     = data['cnm']
	nhi_hl      = data['nhi_hl']
	nh_pl       = data['nh_pl']
	nhi_fk      = data['nhi_fk']

	err_pl      = data['err_pl']
	err_hl      = data['err_hl']
	err_fk      = data['err_fk']
	err_cnm     = cnm_inf['cnm_er']

	diff_fk_hl  = data['fk_hl']
	diff_pl_hl  = data['pl_hl']

	oh          = data['oh']
	src         = data['src']
	idx         = data['idx']

	# Error bar for x-axis and y-axis
	xerr_hl = np.array(err_cnm)/np.array(nhi_cnm)/np.log(10.0) # becaus x-axis = log10(N(HI))
	yerr_pl = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nh_pl, err_pl)

	print yerr_pl

	for k in range(0, len(nhi_cnm)):
		nhi_cnm[k] = np.log10(nhi_cnm[k])

	# Fit and Plot #
	xdata = np.asarray(nhi_cnm)
	ydata = np.asarray(diff_pl_hl)

	# Error bar for x-axis and y-axis
	xerr = np.asarray(xerr_hl)
	yerr = np.asarray(yerr_pl)

	trace = None
	with pm.Model() as model:
	    alpha = pm.Uniform('alpha', lower=0., upper=5.)
	    beta  = pm.Uniform('beta', lower=0., upper=5.)
	    # sigma = pm.Uniform('sigma', lower=0, upper=20)
	    sigma = yerr
	    
	    y_est = alpha + beta * xdata
	    
	    likelihood = pm.Normal('y', mu=y_est, sd=sigma, observed=ydata)
	    
	    # obtain starting values via MAP
	    start = pm.find_MAP()
	    step  = pm.NUTS(state=start)
	    trace = pm.sample(2000, step, start=start, progressbar=False)
	    
	    # pm.traceplot(trace)

	# pprint(trace['alpha'].mean())
	# pprint(trace['alpha'].std())
	# print pm.summary(trace)
	# print pm.summary(trace, ['alpha'])
	# print pm.stats()
	# print(trace.__dict__)

	a = trace['alpha']
	b = trace['beta']

	print a.mean()
	print b.mean()
	xfit = np.linspace(xdata.min(), xdata.max(), 20)
	yfit = b[:, None] * xfit + a[:, None]
	mu   = yfit.mean(0)
	sig  = 1.0*yfit.std(0)

	# MCMC Linear Fit & Plot #
	m  = round(trace['beta'].mean(),2)
	b  = round(trace['alpha'].mean(),2)
	ea = round(trace['beta'].std(),2)
	eb = round(trace['alpha'].std(),2)

	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='b', marker='^', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='Ratio $f = N_{H}$/$N_{HI}$')
	plt.plot(xfit, mu, '-r', mew=2, linewidth=3, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MCMC linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

	plt.title('Correlation between ratio $f = N_{H}$/$N_{HI}$ and Cold Neutral Medium column density $N_{CNM}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	plt.ylabel('$Ratio f = N_{H}$/$N_{HI}$', fontsize=35)
	plt.xlabel('log$_{10}$($N_{CNM}/10^{20}$ cm$^{-2}$)', fontsize=35)
	plt.xlim(-0.6, 1.28)
	plt.ylim(0.5, 3.5)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(-0.5, 0.65, '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(-0.5, 0.55, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(b)+'\pm'+str(eb)+']$', fontsize=20 )

	# plt.text(-0.5, 0.4, 'a = '+str(m)+'$\pm$'+str(ea) +',  b = '+str(b)+'$\pm$'+str(eb), color='blue', fontsize=17)
	# plt.text(-0.5, 0.31, '(Available sources with the presence of OH line are shown)', color='red', fontsize=17)

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=18,
	               )
	plt.show()

## Calculate the uncertainties of factors #
 #
 # params 
 # params 
 # params 
 #
 # return list of factor_uncertainties
 # 
 # Author Van Hiep ##
def uncertainty_of_factors(factr, nhi_hl, err_hl, nhi, err):
	d1 = np.array(err_hl)/np.array(nhi_hl)
	d1 = np.square(d1)

	d2 = np.array(err)/np.array(nhi)
	d2 = np.square(d2)

	d  = np.sqrt(d1+d2)*np.array(factr)

	return d.tolist()

#================= MAIN ========================#
col_density = read_nhi_fukui_nh_planck('../result/26src_no_co_nhi_and_uncertainties_full.txt')
lownhi      = read_23src_lownhi(fname = '../result/23src_lownhi_nhi_and_uncertainties.txt')
plot_planck_factor_vs_nhi(col_density, lownhi) ## Note: just use only one function at once
# plot_planck_factor_vs_nhi_cnm(col_density) ## Note: just use only one function at once 