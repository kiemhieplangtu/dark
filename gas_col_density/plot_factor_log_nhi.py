import os, sys, shutil
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import pymc3             as pm
import copy

from numpy   import array
from restore import restore
import operator

## Read N(HI)_Heiles for each component #
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
	fmt   = ['i',   's',  'f',    'f',      'f',    'f',     'f',    'f',       'f',     'f',   'f',   'f',    'f',       'i',    'f']
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
def read_cnm_err(fname = 'result/26src_no_co_cnm_uncertainties_arcb.txt'):
	cols  = ['idx','src','nhi','nhi_er','cnm','cnm_er','wnm','wnm_er']
	fmt   = ['i',   's',  'f', 'f',     'f',    'f',     'f',     'f']
	dat   = restore(fname, 2, cols, fmt)
	inf   = dat.read()
	return inf

## Plot the correlation between Factors and log(N(HI)) with uncertainties #
 #
 # params dict col_density N(HI) of 26 sources
 #
 # return Void
 # 
 # Author Van Hiep ##
def plot_factor_vs_nhi_err(col_density):
	nhi_hl      = col_density['nhi_hl']
	nh_pl       = col_density['nh_pl']
	nhi_fk      = col_density['nhi_fk']

	err_pl      = col_density['err_pl']
	err_hl      = col_density['err_hl']
	err_fk      = col_density['err_fk']

	diff_fk_hl  = col_density['fk_hl']
	diff_pl_hl  = col_density['pl_hl']

	# Error bar for x-axis and y-axis
	xerr_hl = np.array(err_hl)/np.array(nhi_hl)/np.log(10.0)
	yerr_pl = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nh_pl, err_pl)
	yerr_fk = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nhi_fk, err_fk) # diff_pl_hl just to check, because no error for fk_conversion factor

	for k in range(0, len(nhi_hl)):
		nhi_hl[k] = np.log10(nhi_hl[k])

	# Fit and Plot #
	m, b = np.polyfit(nhi_hl,diff_pl_hl,1)	
	n, c = np.polyfit(nhi_hl,diff_fk_hl,1)

	plt.errorbar(nhi_hl, diff_fk_hl, yerr=yerr_fk, xerr=xerr_hl, fmt='r.', label='')
	plt.errorbar(nhi_hl, diff_pl_hl, yerr=yerr_pl, xerr=xerr_hl, fmt='b.', label='')

	# plt.plot(nhi_hl, diff_fk_hl, 'r.', label='')
	# plt.plot(nhi_hl, diff_pl_hl, 'b.', label='')
	plt.plot(nhi_hl, m*np.array(nhi_hl) + b, 'k-', label='')
	plt.plot(nhi_hl, n*np.array(nhi_hl) + c, 'r-', label='')

	m = round(m, 2)
	b = round(b, 2)
	n = round(n, 2)
	c = round(c, 2)

	plt.xlabel('log(N(HI)_heiles) [1e20/cm2]')
	plt.ylabel('Factor')
	plt.title('Factor vs N(HI)_Heiles')
	plt.grid(True)
	plt.xlim(0, 1.6)
	plt.ylim(0, 5)

	plt.text(0.21, 4.1, 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=12)
	plt.text(0.21, 0.4, 'a = '+str(m)+'  b = '+str(b), color='blue', fontsize=12)
	plt.text(0.21, 2.8, 'a = '+str(n)+'  b = '+str(c), color='red', fontsize=12)

	#plt.legend(loc='upper right')
	plt.show()

## Plot the correlation between Factors and log(N(HI)) #
 #
 # params dict data N(HI) of 26 sources
 #
 # return Void
 # 
 # Author Van Hiep ##
def plot_factor_vs_nhi(data):

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

	# Error bar for x-axis and y-axis
	err_hl  = np.array(err_hl)/np.array(nhi_hl)/np.log(10.0)
	yerr_pl = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nh_pl, err_pl)
	yerr_fk = uncertainty_of_factors(diff_pl_hl, nhi_hl, err_hl, nhi_fk, err_fk)

	for k in range(0, len(nhi_hl)):
		nhi_hl[k] = np.log10(nhi_hl[k])

	# Fit and Plot #
	m, b = np.polyfit(nhi_hl,diff_pl_hl,1)	
	n, c = np.polyfit(nhi_hl,diff_fk_hl,1)

	plt.errorbar(nhi_hl, diff_fk_hl, yerr=yerr_fk, xerr=err_hl, fmt='r.', label='')
	plt.errorbar(nhi_hl, diff_pl_hl, yerr=yerr_pl, xerr=err_hl, fmt='b.', label='')

	# plt.plot(nhi_hl, diff_fk_hl, 'r.', label='')
	# plt.plot(nhi_hl, diff_pl_hl, 'b.', label='')
	plt.plot(nhi_hl, m*np.array(nhi_hl) + b, 'k-', label='')
	plt.plot(nhi_hl, n*np.array(nhi_hl) + c, 'r-', label='')

	m = round(m, 2)
	b = round(b, 2)
	n = round(n, 2)
	c = round(c, 2)

	plt.xlabel('log(N(HI)_heiles) [1e20/cm2]')
	plt.ylabel('Factor')
	plt.title('Factor vs N(HI)_Heiles')
	plt.grid(True)
	plt.xlim(0,1.6)
	plt.ylim(0,5)

	plt.text(0.21, 4.1, 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=12)
	plt.text(0.21, 0.4, 'a = '+str(m)+'  b = '+str(b), color='blue', fontsize=12)
	plt.text(0.21, 2.8, 'a = '+str(n)+'  b = '+str(c), color='red', fontsize=12)

	for i,j in zip(nhi_hl,diff_pl_hl):
		if (j>-1.0) :
			plt.annotate('('+str(j)+')',xy=(i,j))

	#plt.legend(loc='upper right')
	plt.show()

## Plot the correlation between Factors and log(N(HI)) only from Planck #
 #
 # params dict data N(HI) of 26 sources
 #
 # return Void
 # 
 # Author Van Hiep ##
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
	xdata = np.asarray(nhi_hl)
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
	    step = pm.NUTS(state=start)
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

	plt.title('Correlation between $f = N_{H}$/$N_{HI}$ and HI column density $N_{HI}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	plt.ylabel('$Ratio f = N_{H}$/$N_{HI}$', fontsize=35)
	plt.xlabel('log$_{10}$($N_{HI}/10^{20}$ cm$^{-2}$)', fontsize=35)
	plt.xlim(0, 1.6)
	plt.ylim(0, 3)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(0.2, 0.31, '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(0.2, 0.4, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(b)+'\pm'+str(eb)+']$', fontsize=20 )

	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(src)):
		if (oh[i] > 0) :
			plt.annotate('('+str(src[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
	               xytext=(-50.,30.), textcoords='offset points',
	               arrowprops=dict(arrowstyle="->"),fontsize=18,
	               )
	plt.show()

## Plot the correlation between Factors and log(N(HI)_CNM) only from Planck #
 #
 # params dict data N(HI) of 26 sources
 #
 # return Void
 # 
 # Author Van Hiep ##
def plot_planck_factor_vs_nhi_cnm(data):

	cnm_inf     = read_cnm_err() # the Uncertainties of CNM component

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

	dif_pl_cnm = list_div(nh_pl, nhi_cnm)

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

	plt.title('Correlation between Factor $f = N_{H}$/$N_{HI}$ and Cold Neutral Medium column density $N_{CNM}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	plt.ylabel('$Ratio f = N_{H}$/$N_{HI}$', fontsize=35)
	plt.xlabel('log$_{10}$($N_{CNM}/10^{20}$ cm$^{-2}$)', fontsize=35)
	plt.xlim(-0.6, 1.28)
	plt.ylim(0., 3.)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	plt.text(-0.5, 0.31, '(Available sources with the presence of OH are shown)', color='k', fontsize=17)
	plt.text(-0.5, 0.4, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(b)+'\pm'+str(eb)+']$', fontsize=20 )

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

	# # Fit and Plot #

	# plt.errorbar(nhi_cnm, diff_pl_hl, yerr=yerr_pl, xerr=xerr_hl, fmt='ro', label='')
	# # plt.plot(nhi_cnm, diff_pl_hl, 'ro', label='Factor $f = N_{H}$/$N_{HI}$')
	# plt.plot(nhi_cnm, diff_pl_hl, 'ro', label='Factor $f = N_{H}$/$N_{HI}$')
	# plt.plot(nhi_cnm, m*np.array(nhi_cnm) + b, 'b-', label='Best linear fit', linewidth=2)

	# m  = round(m,2)
	# b  = round(b,2)
	# ea = round(ea,2)
	# eb = round(eb,2)

	# plt.ylabel('$Factor f = N_{H}$/$N_{HI}$', fontsize=35)
	# plt.xlabel('log$_{10}$($N_{CNM}/10^{20}$ cm$^{-2}$)', fontsize=35)
	# plt.title('Correlation between Factor $f = N_{H}$/$N_{HI}$ and CNM column density $N_{CNM}$ \nalong 26 lines-of-sight without the presence of CO line', fontsize=30)
	# plt.grid(True)
	# plt.tick_params(axis='x', labelsize=18)
	# plt.tick_params(axis='y', labelsize=18)

	

	

	# for i,j in zip(idx,diff_pl_hl):
	# 	j = j + 0.02
	# 	if (oh[i] > 0) :
	# 		plt.annotate('('+str(src[i])+')',xy=(nhi_cnm[i],j), fontsize=12)

	# plt.legend(loc='upper left', fontsize=18)
	# plt.show()

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

## plot N(HI)/N(HI)_thin vs N(HI)_thin #
 #
 # params dict data data to plot
 #
 # return void
 # 
 # Author Van Hiep ##
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

## linear fit #
 #
 # params x list x-data
 # params y list y-data
 #
 # return fit parameters and errors
 # 
 # Author Van Hiep ##
def linear_fit(x,y):
	ret ={}

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

	ret['a']=a
	ret['b']=b
	ret['ea']=er_a
	ret['eb']=er_b

	return ret

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

col_density = read_nhi_fukui_nh_planck()
# plot_planck_factor_vs_nhi(col_density) ## Note: just use only one function at once
plot_planck_factor_vs_nhi_cnm(col_density) ## Note: just use only one function at once