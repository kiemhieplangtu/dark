import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import pymc3             as pm
import operator

from numpy        import array
from restore      import restore
from plotting     import cplot
from mpfit        import mpfit
import copy

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

## plot N(OH) vs N(HI) #
 # to check if N(H2) = 1e7*N(OH)
 #
 # params None
 # return void
 # 
 # version 09/2016 
 # author Nguyen Van Hiep ##
def plot_oh_vs_h2():
	## N(H) N(HI) N(OH) from 19 sources ##
	cols = ['idx','src', 'noh','noh_er', 'nhi','nhi_er', 'nh','nh_er']
	fmt  = ['i',  's',    'f',  'f',      'f',   'f',     'f',  'f']
	data = restore('sub_data/nh_nhi_noh_19src_er.txt', 3, cols, fmt)
	dat  = data.read()

	sc     = dat['src']
	idx    = dat['idx']
	nhi    = dat['nhi']
	nhi_er = dat['nhi_er']
	nh     = dat['nh']
	nh_er  = dat['nh_er']
	noh    = dat['noh']
	noh_er = dat['noh_er']

	nhi    = np.asarray(dat['nhi'])
	nhi_er = np.asarray(dat['nhi_er'])
	nh     = np.asarray(dat['nh'])
	nh_er  = np.asarray(dat['nh_er'])
	noh    = np.asarray(dat['noh'])
	noh_er = np.asarray(dat['noh_er'])

	xdata = noh*1e14
	ydata = nh*1e20

	# Error bar for x-axis and y-axis
	xerr = noh_er*1e14
	yerr = nh_er*1e20

	## Fit  Tau, V0, Width ##
	lguess  = [ 2.5e7, 1.e21]

	npar    =  len(lguess)
	guessp  = np.array(lguess, dtype='float64')
	plimd   = [[False,False]]*npar
	plims   = [[0.,0.]]*npar
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
	pname   = ['slope','offset']
	pfix    = [False]*npar

	parinfo = []
	for i in range(len(guessp)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		parinfo[i]['value']   = guessp[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]

	##  1665 ###
	x    = xdata.astype(np.float64)
	y    = ydata.astype(np.float64)
	er   = yerr.astype(np.float64)

	fa = {'x':x, 'y':y, 'err':er}
	mp = mpfit(myfunc, guessp, parinfo=parinfo, functkw=fa, quiet=True)

	## ********* Results ********* ##
	print '********* Results *********'
	abp   = mp.params
	abper = mp.perror
	for i in range(len(parinfo)):
		print "%s = %f +/- %f" % (parinfo[i]['parname'],abp[i],abper[i])
	## Plot ##
	a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
	b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
	xfit  = np.linspace(xdata.min(), xdata.max(), 20)
	yfit  = a[:, None] * xfit + b[:, None]
	mu    = yfit.mean(0)
	sig   = 1.0*yfit.std(0)
	fit   = abp[0]*x+abp[1]

	m  = round(abp[0]/1e7,2)
	b  = round(abp[1]/1e21,2)
	ea = round(abper[0]/1e7,2)
	eb = round(abper[1]/1e21,2)

	plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$N(H_{2}) = N(H)_{Dust} - N(HI)_{MS}$')
	plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
	# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
	# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

	plt.title('$N(H_{2}) vs N(OH)$', fontsize=30)
	plt.ylabel('$N(H_{2})$', fontsize=35)
	plt.xlabel('$N(OH)$', fontsize=35)
	plt.xlim(0.,5.5e14)
	plt.ylim(-2e21,1.9e22)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	plt.text(5e13,-4e20, 'a = ['+str(m)+'$\pm$'+str(ea) +']$10^{7}$,  b = ['+str(b)+'$\pm$'+str(eb)+']$10^{21}$', color='blue', fontsize=17)
	plt.text(2.8e14,0.0, 'N(H) from dust: $tau_{353} = [8.4\pm3.0]\cdot 10^{-27}N_{H} $', color='blue', fontsize=17)
	plt.text(2.8e14,0.05e22, 'Previous results: $N_{H_{2}} = 10^{7}N_{OH} $', color='blue', fontsize=17)
	plt.text(2.8e14,0.1e22, 'My result: $N_{H_{2}} = '+str(m)+'\cdot 10^{7}N_{OH}+'+str(b)+'\cdot 10^{21}$', color='blue', fontsize=17)
	plt.legend(loc='lower right', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	for i in range(len(sc)):
		plt.annotate('('+str(sc[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
                xytext=(-50.,30.), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),fontsize=18,
                )
	plt.show()

#================= MAIN ========================#
plot_oh_vs_h2()