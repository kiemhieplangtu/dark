import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
import copy

from   restore           import restore
from   mpfit             import mpfit

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
		sx  = sx  + x[i]
		sy  = sy  + y[i]
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

	ret['a']  = a
	ret['b']  = b
	ret['ea'] = er_a
	ret['eb'] = er_b

	return ret

## Read info of 94 SPONGE + MS sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_info_94src(fname = '../rearrange/nhi_lb_thin_78src.txt'):
	cols = ['idx','src','l', 'b', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_info_no_co(fname = '../../co12/result/26src_no_co_with_sponge.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j', 'oh', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',    'i', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Read info of each source to dictionary #
 #
 # params dict dat Data
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read2dict(dat):
	sc     = dat['src']
	xl     = dat['l']
	xb     = dat['b']
	hi     = dat['nhi']
	hier   = dat['nhi_er']
	thin   = dat['thin']
	thiner = dat['thin_er']
	cnm    = dat['cnm']
	cnmer  = dat['cnm_er']
	wnm    = dat['wnm']
	wnmer  = dat['wnm_er']

	ret = {}
	for i in range(len(sc)):
		ret[sc[i]]           = {}
		ret[sc[i]]['l']      = xl[i]
		ret[sc[i]]['b']      = xb[i]
		ret[sc[i]]['nhi']    = hi[i]
		ret[sc[i]]['nhi_er'] = hier[i]
		ret[sc[i]]['thin']   = thin[i]
		ret[sc[i]]['thiner'] = thiner[i]
		ret[sc[i]]['cnm']    = cnm[i]
		ret[sc[i]]['cnm_er'] = cnmer[i]
		ret[sc[i]]['wnm']    = wnm[i]
		ret[sc[i]]['wnm_er'] = wnmer[i]

	return ret

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

## ======== MAIN ============= ##

## 94 MS+SPONGE Sources, MS prior
sc94   = read_info_94src(fname = '../result/nhi_thin_cnm_wnm_94src.txt')
# sc     = sc94['src']
# xl     = sc94['l']
# xb     = sc94['b']
# hi     = sc94['nhi']
# hier   = sc94['nhi_er']
# thin   = sc94['thin']
# thiner = sc94['thin_er']
# cnm    = sc94['cnm']
# cnmer  = sc94['cnm_er']
# wnm    = sc94['wnm']
# wnmer  = sc94['wnm_er']

## Consider |b| > 10 Only
sc     = []
xl     = []
xb     = []
hi     = []
hier   = []
thin   = []
thiner = []
cnm    = []
cnmer  = []
wnm    = []
wnmer  = []
for i in range(len(sc94['src'])):
	if(np.abs(sc94['b'][i]) > 10.):
		sc.append(sc94['src'][i])
		xl.append(sc94['l'][i])
		xb.append(sc94['b'][i])
		hi.append(sc94['nhi'][i])
		hier.append(sc94['nhi_er'][i])
		thin.append(sc94['thin'][i])
		thiner.append(sc94['thin_er'][i])
		cnm.append(sc94['cnm'][i])
		cnmer.append(sc94['cnm_er'][i])
		wnm.append(sc94['wnm'][i])
		wnmer.append(sc94['wnm_er'][i])
		
# Fit and Plot #
params = linear_fit(thin,hi)
a      = params['a']
b      = params['b']
ea     = params['ea']
eb     = params['eb']

# plt.plot(thin, hi, 'b^', ms=10, label='$N_{HI} vs N^{thin}_{HI}$')
plt.errorbar(thin, hi, yerr=hier, xerr=thiner, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$N_{HI} vs N^{thin}_{HI}$')
plt.plot(thin, a*np.array(thin) + b, 'r-', linewidth=4, label='Best linear fit')

a  = round(a, 2)
b  = round(b, 2)
ea = round(ea, 2)
eb = round(eb, 2)

plt.grid()
plt.title('$N_{HI}$ vs $N^{thin}_{HI}$, 21-SPONGE + MS*', fontsize = 35)
plt.ylabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize = 35)
plt.xlabel('$N^{thin}_{HI} [10^{20} cm^{-2}]$', fontsize = 35)
plt.xlim(-1., 50.)
plt.ylim(-1., 170.)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(10., 5., '$N_{HI}[10^{20} cm^{-2}] = ('+str(a)+'\pm'+str(ea) +')\cdot N^*_{HI}[10^{20} cm^{-2}] + ('+str(b)+'\pm'+str(eb)+')$', color='blue', fontsize=20)
plt.legend(loc='upper right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')

plt.show()




## 94 MS+SPONGE Sources, 21SPONGE* prior
sc94 = read_info_94src(fname = '../result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
# sc     = sc94['src']
# xl     = sc94['l']
# xb     = sc94['b']
# hi     = sc94['nhi']
# hier   = sc94['nhi_er']
# thin   = sc94['thin']
# thiner = sc94['thin_er']
# cnm    = sc94['cnm']
# cnmer  = sc94['cnm_er']
# wnm    = sc94['wnm']
# wnmer  = sc94['wnm_er']
sc     = []
xl     = []
xb     = []
hi     = []
hier   = []
thin   = []
thiner = []
cnm    = []
cnmer  = []
wnm    = []
wnmer  = []
for i in range(len(sc94['src'])):
	if(np.abs(sc94['b'][i]) > 10.):
		sc.append(sc94['src'][i])
		xl.append(sc94['l'][i])
		xb.append(sc94['b'][i])
		hi.append(sc94['nhi'][i])
		hier.append(sc94['nhi_er'][i])
		thin.append(sc94['thin'][i])
		thiner.append(sc94['thin_er'][i])

# Fit and Plot #
params = linear_fit(thin,hi)
a      = params['a']
b      = params['b']
ea     = params['ea']
eb     = params['eb']

# plt.plot(thin, hi, 'b^', ms=10, label='$N_{HI} vs N^{thin}_{HI}$')
plt.errorbar(thin, hi, yerr=hier, xerr=thiner, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$N_{HI} vs N^{thin}_{HI}$')
plt.plot(thin, a*np.array(thin) + b, 'r-', linewidth=4, label='Best linear fit')

a  = round(a, 2)
b  = round(b, 2)
ea = round(ea, 2)
eb = round(eb, 2)

plt.grid()
plt.title('$N_{HI}$ vs $N^{thin}_{HI}$, 21-SPONGE* + MS', fontsize = 35)
plt.ylabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize = 35)
plt.xlabel('$N^{thin}_{HI} [10^{20} cm^{-2}]$', fontsize = 35)
plt.xlim(-1., 50.)
plt.ylim(-1., 170.)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(10., 5., '$N_{HI}[10^{20} cm^{-2}] = ('+str(a)+'\pm'+str(ea) +')\cdot N^*_{HI}[10^{20} cm^{-2}] + ('+str(b)+'\pm'+str(eb)+')$', color='blue', fontsize=20)
plt.legend(loc='upper right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')

plt.show()


## 94 MS+SPONGE Sources, 21SPONGE* prior, R vs log10(N*HI)
fact   = []
lognhi = []
for i in range(0, len(thin)):
	temp = round(hi[i]/thin[i], 3)
	lognhi.append(np.log10(thin[i]))
	fact.append(temp)

# Error bar for x-axis and y-axis
thin   = np.asarray(thin)
hi     = np.asarray(hi)
xdata  = np.asarray(lognhi)
ydata  = np.asarray(fact)
hier   = np.array(hier)
thiner = np.array(thiner)
xerr   = thiner/thin/np.log(10.0)
yerr   = uncertainty_of_factors(ydata, thin, thiner, hi, hier)

# Fit and Plot #
params = linear_fit(lognhi,fact)
a      = params['a']
b      = params['b']
ea     = params['ea']
eb     = params['eb']

# plt.plot(lognhi, fact, 'b^', label='Ratio $f = N_{HI}$/$N^*_{HI}$', markersize=10)
plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='b', marker='o', ls='None', markersize=4, markeredgecolor='b', markeredgewidth=1, label='Ratio $f = N_{HI}$/$N^*_{HI}$')
plt.plot(lognhi, a*np.array(lognhi) + b, 'r-', linewidth=4, label='Best linear fit')

a  = round(a, 2)
b  = round(b, 2)
ea = round(ea, 2)
eb = round(eb, 2)

plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
plt.title('Correlation between Total HI column densities $N_{HI}$ and \n HI optically thin column densities $N^*_{HI}$ along 94 (Millennium Survey & 21-SPONGE) lines-of-sight', fontsize=30)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)

plt.text(0.0, 3.0, '$f = ['+str(a)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)
plt.text(0.0, 3.2, r'$f = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al.', color='blue', fontsize=20)

plt.legend(loc='upper left', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')
plt.show()

########### MPFIT ############
xdata = xdata
ydata = ydata

# Error bar for x-axis and y-axis
xerr = xerr
yerr = yerr

## Fit ##
lguess  = [ 0.33, 0.85]

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

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata,ydata,xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='Ratio $f = N_{HI}$/$N^*_{HI}$')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
# plt.plot(xfit, alpha*xfit[:, None] + beta, c='blue', alpha=0.01)
# plt.plot(xfit, alpha_stats['mean']*xfit + beta_stats['mean'], linewidth=2, c='red')

plt.title('Correlation between Total HI column densities $N_{HI}$ and \n HI optically thin column densities $N^*_{HI}$ along 94 (Millennium Survey & 21-SPONGE) lines-of-sight', fontsize=30)
plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
plt.xlim(0.0, 2.0)
plt.ylim(-1.0, 6.0)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)

plt.text(0.0, 3.0, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)
plt.text(0.0, 3.2, r'$f = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al.', color='blue', fontsize=20)
plt.legend(loc='lower right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')
# for i in range(len(sc)):
# 	plt.annotate('('+str(sc[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
#             xytext=(-50.,30.), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"),fontsize=18,
#             )
plt.show()
########### END - MPFIT ############




## Check Low NHI sources, NHI < 3.0 e20 cm-2
## 94 MS+SPONGE Sources, 21SPONGE* prior
sc94 = read_info_94src(fname = '../result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
sc     = sc94['src']
xl     = sc94['l']
xb     = sc94['b']
hi     = sc94['nhi']
hier   = sc94['nhi_er']
thin   = sc94['thin']
thiner = sc94['thin_er']
cnm    = sc94['cnm']
cnmer  = sc94['cnm_er']
wnm    = sc94['wnm']
wnmer  = sc94['wnm_er']
src26 = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
sc26  = src26['src']
k     = 0
for i in range(len(sc)):
	if ( (hi[i] < 3.) and (sc[i] not in sc26) ):
		k = k + 1
		print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
			.format(k-1, sc[i], xl[i], xb[i], hi[i], hier[i], thin[i], thiner[i], cnm[i], cnmer[i], wnm[i], wnmer[i]  ))         ## lownhi_thin_cnm_wnm.txt