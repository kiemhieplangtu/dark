import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
from   restore           import restore

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
def read_info_94src(fname = '../result/nhi_thin_cnm_wnm_94src.txt'):
	cols = ['idx','src','l', 'b', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict infocd 
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

## ======== MAIN ============= ##

## 94 MS+SPONGE Sources
sc94 = read_info_94src(fname = '../result/nhi_thin_cnm_wnm_94src.txt')
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
plt.title('$N_{HI}$ vs $N^{thin}_{HI}$, 21SPONGE + MS', fontsize = 35)
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

## Check Low NHI sources, NHI < 3.0 e20 cm-2
info  = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
sc26  = src26['src']
k     = 0
for i in range(len(sc)):
	if ( (hi[i] < 3.) and (sc[i] not in sc26) ):
		k = k + 1
		print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
			.format(k-1, sc[i], xl[i], xb[i], hi[i], hier[i], thin[i], thiner[i], cnm[i], cnmer[i], wnm[i], wnmer[i]  ))         ## lownhi_thin_cnm_wnm.txt