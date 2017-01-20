import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy        import array
from restore      import restore
from plotting     import cplot
import pymc3 as pm
import copy
from mpfit import mpfit

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
	## N(H) N(HI) N(OH) from 19 sources, created by cal_noh_nh_nhi.py ##
	cols = ['idx','src', 'noh','noh_er', 'nhi','nhi_er', 'nh','nh_er']
	fmt  = ['i',  's',    'f',  'f',      'f',   'f',     'f',  'f']
	data = restore('sub_data/nh_nhi_noh_19src_er.txt', 3, cols, fmt)
	dat  = data.read()

	src    = dat['src']
	idx    = dat['idx']
	nhi    = dat['nhi']
	nhi_er = dat['nhi_er']
	nh     = dat['nh']
	nh_er  = dat['nh_er']
	noh    = dat['noh']
	noh_er = dat['noh_er']

	xnoh  = []
	hdiff = [] # hdiff = N(H2) = N(H) - N(HI)
	h2_er = []
	x_er  = []
	sc    = []
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

	# Error bar for x-axis and y-axis
	xerr = x_er*1e14
	yerr = h2_er*1e20

	## Fit ##
	lguess  = [ 1.e7, 1.e21]

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

	# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label=r'$N(H_{2}) = \frac{1}{2}[N(H)_{Dust} - N(HI)_{MS}]$')
	plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
	# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
	# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

	plt.title('$N(H_{2}) vs N(OH)$', fontsize=30)
	plt.ylabel('$N(H_{2})$', fontsize=35)
	plt.xlabel('$N(OH)$', fontsize=35)
	plt.xlim(0.0,4.0e14)
	plt.ylim(-1.0e21,5.0e21)
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

## plot N(H) from Dust vs N(HI) for 19 src with OH#
 # to check if N(H2) = 1e7*N(OH)
 #
 # params None
 # return void
 # 
 # version 10/2016 
 # author Nguyen Van Hiep ##
def plot_nh_vs_nhi():
	## N(H) N(HI) N(OH) from 19 sources ##
	cols = ['idx','src', 'noh','noh_er', 'nhi','nhi_er', 'nh','nh_er']
	fmt  = ['i',  's',    'f',  'f',      'f',   'f',     'f',  'f']
	data = restore('sub_data/nh_nhi_noh_19src_er.txt', 3, cols, fmt)
	dat  = data.read()

	src    = dat['src']
	idx    = dat['idx']
	nhi    = dat['nhi']
	nhi_er = dat['nhi_er']
	nh     = dat['nh']
	nh_er  = dat['nh_er']
	noh    = dat['noh']
	noh_er = dat['noh_er']

	# ii = src.index('T0629+10') # high N(H) l = 201.53, b=0.5079 in Gal Plane

	# nhi.pop(ii)
	# nhi_er.pop(ii)
	# nh.pop(ii)
	# nh_er.pop(ii)
	# noh.pop(ii)
	# noh_er.pop(ii)

	nhi    = np.asarray(nhi)
	nhi_er = np.asarray(nhi_er)
	nh     = np.asarray(nh)
	nh_er  = np.asarray(nh_er)
	noh    = np.asarray(noh)
	noh_er = np.asarray(noh_er)
	nh2    = noh*10. ## N(H2) = N(OH)*1e7

	xdata = np.log10(nhi)
	ydata = nh/nhi

	print ydata

	# Error bar for x-axis and y-axis
	xerr = nhi_er/nhi/np.log(10.0)
	yerr = uncertainty_of_factors(ydata, nhi, nhi_er, nh, nh_er)

	## Fit ##
	lguess  = [ 0.3, 1.0]

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

	m  = round(abp[0],2)
	b  = round(abp[1],2)
	ea = round(abper[0],2)
	eb = round(abper[1],2)

	plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$N(H_{2}) = N(H)_{Dust} - N(HI)_{MS}$')
	plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
	# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
	# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

	plt.title('$N(H_{2}) vs N(OH)$', fontsize=30)
	plt.ylabel('$N(H_{2})$', fontsize=35)
	plt.xlabel('$N(OH)$', fontsize=35)
	plt.xlim(0.,2.5)
	# plt.ylim(0.,11.0)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	plt.text(0.01,9.0, '$tau_{353} = [8.4\pm3.0]10^{-27}N_{H} $', color='k', fontsize=15)
	plt.text(0.01,8.5, '$N_{H_{2}} = [1.0\pm??]10^{7}N_{OH} $', color='k', fontsize=15)
	plt.text(0.01,4.0, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(b)+'\pm'+str(eb)+']\cdot$', color='k', fontsize=20)
	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(sc)):
	# 	plt.annotate('('+str(sc[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
 #                xytext=(-50.,30.), textcoords='offset points',
 #                arrowprops=dict(arrowstyle="->"),fontsize=18,
 #                )
	plt.show()

## plot ([N(H)_from_Dust] vs [N(HI) + N(H2)_from_OH] ) for 19 src with OH#
 # to check if N(H2) = 1e7*N(OH)
 #
 # params None
 # return void
 # 
 # version 10/2016 
 # author Nguyen Van Hiep ##
def plot_estimate_nh_vs_nhi_nh2():
	## N(H) N(HI) N(OH) from 19 sources ##
	cols = ['idx','src', 'noh','noh_er', 'nhi','nhi_er', 'nh','nh_er']
	fmt  = ['i',  's',    'f',  'f',      'f',   'f',     'f',  'f']
	data = restore('sub_data/nh_nhi_noh_19src_er.txt', 3, cols, fmt)
	dat  = data.read()

	src    = dat['src']
	idx    = dat['idx']
	nhi    = dat['nhi']
	nhi_er = dat['nhi_er']
	nh     = dat['nh']
	nh_er  = dat['nh_er']
	noh    = dat['noh']
	noh_er = dat['noh_er']

	ii = src.index('T0629+10') # high N(H) l = 201.53, b=0.5079 in Gal Plane

	nhi.pop(ii)
	nhi_er.pop(ii)
	nh.pop(ii)
	nh_er.pop(ii)
	noh.pop(ii)
	noh_er.pop(ii)

	nhi    = np.asarray(nhi)
	nhi_er = np.asarray(nhi_er)
	nh     = np.asarray(nh)
	nh_er  = np.asarray(nh_er)
	noh    = np.asarray(noh)
	noh_er = np.asarray(noh_er)
	nh2    = noh*10. ## N(H2) = N(OH)*1e7

	# for i in range(len(src)):
	# 	print ("{}   {}\t{}\t{}\t{}\t{}"\
	# 	 .format(i, src[i],noh[i],nhi[i],nh[i],round(nh2[i],3)   ) )

	xdata = np.log10(nhi+nh2)
	ydata = nh/(nhi+nh2)

	# Error bar for x-axis and y-axis
	xerr  = nhi_er/(nhi+nh2)/np.log(10.0)
	yerr  = uncertainty_of_factors(ydata, nhi, nhi_er, nh, nh_er)

	## MCMC pymc3 linear fit & Plot, y = alpha + beta * x ##
	lfit = pymc3lfit()
	a,b  = lfit.fit(xdata,ydata,yerr,arange=[0.,5.0],brange=[0.,2.0] )

	print a.mean()
	print b.mean()
	xfit = np.linspace(xdata.min(), xdata.max(), 20)
	yfit = b[:, None] * xfit + a[:, None]
	mu   = yfit.mean(0)
	sig  = 1.0*yfit.std(0)

	beta  = round(b.mean(),2)
	alpha = round(a.mean(),2)
	berr  = round(b.std(),2)
	aerr  = round(a.std(),2)

	# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Ratio $f = N_{HI}$/$N^*_{HI}$')
	plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='Ratio $f = N_{H}/(N_{HI}+N_{H_{2}})$, $N_{H}$ from thermal dust, $N_{H_{2}}$ from OH')
	plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='Linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
	# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
	# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

	plt.title('Correlation between Ratio $f = N_{H}$/$(N_{HI}+N_{H_{2}})$ and HI column density $(N_{HI}+N_{H_{2}})$ \nalong 19 lines-of-sight with OH', fontsize=30)
	plt.ylabel('$Ratio f = N_{H}$/$(N_{HI}+N_{H_{2}})$', fontsize=35)
	plt.xlabel('$log_{10}((N_{HI}+N_{H_{2}})/10^{20} cm^{-2}$)', fontsize=35)
	plt.xlim(0.,2.5)
	plt.ylim(0.,4.)
	plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)

	plt.text(0.01,9.0, '$tau_{353} = [8.4\pm3.0]10^{-27}N_{H} $', color='k', fontsize=15)
	plt.text(0.01,8.5, '$N_{H_{2}} = [1.0\pm??]10^{7}N_{OH} $', color='k', fontsize=15)
	plt.text(0.01,4.0, '$f = ['+str(beta)+'\pm'+str(berr) +']\cdot log_{10}(N_{HI}/10^{20})+['+str(alpha)+'\pm'+str(aerr)+']\cdot$', color='k', fontsize=20)
	plt.legend(loc='upper left', fontsize=18)
	# plt.savefig("test.png",bbox_inches='tight')
	# for i in range(len(sc)):
	# 	plt.annotate('('+str(sc[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
    #             xytext=(-50.,30.), textcoords='offset points',
    #             arrowprops=dict(arrowstyle="->"),fontsize=18,
    #             )
	plt.show()

#================= MAIN ========================#
plot_oh_vs_h2()
# plot_nh_vs_nhi()
# plot_estimate_nh_vs_nhi_nh2()