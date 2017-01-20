import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import pymc3             as pm
import operator

from numpy           import array
from restore         import restore
from plotting        import cplot
from pymc3linear     import pymc3lfit


## Calculate the uncertainties of factors #
 #
 # params 1D-array factr y-axis
 # params 1D-array nhi N(HI)
 # params 1D-array nhi_er Uncertainties of N(HI)
 # params 1D-array nh N(H)
 # params 1D-array nh_er Uncertainties of N(H)
 #
 # return factor_uncertainties
 # 
 # Author Van Hiep ##
def uncertainty_of_factors(factr, nhi, nhi_er, nh, nh_er):
	d1 = nhi_er/nhi
	d1 = d1**2

	d2 = nh_er/nh
	d2 = d2**2

	d  = np.sqrt(d1+d2)*factr

	return d

## plot N(OH) vs N(HI) #
 # to check if N(H2) = 1e7*N(OH)
 #
 # params None
 # return void
 # 
 # version 09/2016 
 # author Nguyen Van Hiep ##
def plot_oh_vs_h2():
	## N(H) N(HI) N(OH) from 19 sources, created by cal_noh_nh_nhi.py ##
	cols   = ['idx','src', 'noh','noh_er', 'nhi','nhi_er', 'nh','nh_er']
	fmt    = ['i',  's',    'f',  'f',      'f',   'f',     'f',  'f']
	data   = restore('result/ebvnh_nhi_noh_19src_er.txt', 3, cols, fmt)
	dat    = data.read()

	src    = dat['src']
	idx    = dat['idx']
	nhi    = dat['nhi']
	nhi_er = dat['nhi_er']
	nh     = dat['nh']
	nh_er  = dat['nh_er']
	noh    = dat['noh']
	noh_er = dat['noh_er']

	xnoh   = []
	hdiff  = [] # hdiff = N(H2) = N(H) - N(HI)
	h2_er  = []
	x_er   = []
	sc     = []
	for i in range(len(src)):
		if( ((nh[i]-nhi[i])>0.) and ((nh[i]-nhi[i])<200.) ):
		# if( (nh[i]-nhi[i])<200. ):
			xnoh.append(noh[i])
			x_er.append(noh_er[i])
			hdiff.append( (nh[i]-nhi[i])/2. )
			h2_er.append( 0.5*np.sqrt( nh_er[i]**2 + nhi_er[i]**2 ) )
			sc.append(src[i])


	xnoh  = np.asarray(xnoh)
	hdiff = np.asarray(hdiff)
	h2_er = np.asarray(h2_er)
	x_er  = np.asarray(x_er)

	xdata = xnoh*1e14
	ydata = hdiff*1e20 # N(H2)

	aa    = ydata/xdata
	print aa
	print 'Point-to-point N(H)/N(OH)', aa.mean()/1e7

	# Error bar for x-axis and y-axis
	xerr = x_er*1e14
	yerr = h2_er*1e20
	
	## MCMC pymc3 linear fit & Plot, y = alpha + beta * x ##
	lfit = pymc3lfit()
	a,b  = lfit.fit(xdata,ydata,yerr,arange=[0e21,2.0e21],brange=[0e7,2.0e7] )

	print a.mean()
	print b.mean()
	xfit = np.linspace(xdata.min(), xdata.max(), 20)
	yfit = b[:, None] * xfit + a[:, None]
	mu   = yfit.mean(0)
	sig  = 1.0*yfit.std(0)

	beta  = round(b.mean()/1e7,2)
	alpha = round(a.mean()/1e21,2)
	berr  = round(b.std()/1e7,2)
	aerr  = round(a.std()/1e21,2)

	# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Ratio $f = N_{HI}$/$N^*_{HI}$')
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label=r'$N(H_{2}) = \frac{1}{2}[N(H)_{Dust} - N(HI)_{MS}]$')
	plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MCMC linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
	# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
	# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

	plt.title('Correlation between $H_{2}$ column density ($N_{H_{2}}$) and OH column density ($N_{OH}$)', fontsize=30)
	plt.ylabel('$N(H_{2}) [cm^{-2}]$', fontsize=35)
	plt.xlabel('$N(OH) [cm^{-2}]$', fontsize=35)
	plt.xlim(0.0,4.0e14)
	plt.ylim(-1.0e21,5.0e21)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	# plt.text(5e13,-4e20, 'a = ['+str(m)+'$\pm$'+str(ea) +']$10^{7}$,  b = ['+str(b)+'$\pm$'+str(eb)+']$10^{21}$', color='blue', fontsize=17)
	plt.text(2.5e14,0.001e22, '$N_{H}$ obtained from thermal dust: $tau_{353} = [8.4\pm3.0]\cdot 10^{-27}N_{H}$', color='k', fontsize=15)
	# plt.text(2.8e14,0.05e22, 'Previous results: $N_{H_{2}} = 10^{7}N_{OH} $', color='blue', fontsize=17)
	plt.text(2.5e14,0.05e22, '$N_{H_{2}} = ['+str(beta)+'\pm'+str(berr) +']\cdot 10^{7}N_{OH}+['+str(alpha)+'\pm'+str(aerr)+']\cdot 10^{21}$', color='k', fontsize=20)
	# plt.text(2.8e14,0.1e22, 'Linear fit: $N_{H_{2}} = '+str(m)+'\cdot 10^{7}N_{OH}+'+str(b)+'\cdot 10^{21}$', color='blue', fontsize=17)
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