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

## Read info of 30 SPONGE sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt'):
	cols = ['src','l', 'b', 'cnm','cnm_er','wnm','wnm_er','nhi','nhi_er','thin','thin_er']
	fmt  = ['s',  'f', 'f', 'f',    'f',   'f',  'f',     'f',   'f',     'f',    'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Read info of 30 SPONGE sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_info_ms_78src(fname = '../rearrange/nhi_lb_thin_78src.txt'):
	cols = ['idx','src','l', 'b', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
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
## 21SPONGE Claire 30 Sources
sp30sc  = read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt')
sc1     = sp30sc['src']
xl1     = sp30sc['l']
xb1     = sp30sc['b']
hi1     = sp30sc['nhi']
hi1er   = sp30sc['nhi_er']
thin1   = sp30sc['thin']
thin1er = sp30sc['thin_er']
cnm1    = sp30sc['cnm']
cnm1er  = sp30sc['cnm_er']
wnm1    = sp30sc['wnm']
wnm1er  = sp30sc['wnm_er']
spdat   = read2dict(sp30sc)

## 78 MS Sources
ms78sc  = read_info_ms_78src(fname = '../result/nhi_lb_thin_cnm_wnm_78src.txt')
sc2     = ms78sc['src']
xl2     = ms78sc['l']
xb2     = ms78sc['b']
hi2     = ms78sc['nhi']
hi2er   = ms78sc['nhi_er']
thin2   = ms78sc['thin']
thin2er = ms78sc['thin_er']
cnm2    = ms78sc['cnm']
cnm2er  = ms78sc['cnm_er']
wnm2    = ms78sc['wnm']
wnm2er  = ms78sc['wnm_er']
msdat   = read2dict(ms78sc)

## Common & different Sources
comsc = []  ## 14 src
difsc = []  ## 16 src
for sc in spdat:
	if (sc in msdat):
		comsc.append(sc)
	else:
		difsc.append(sc)

for i in range(len(comsc)):
	sc = comsc[i]
	# print sc, spdat[sc]['l'], msdat[sc]['l'], spdat[sc]['l'] - msdat[sc]['l']
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(spdat[sc]['nhi'] - msdat[sc]['nhi'])/msdat[sc]['nhi']

# Fit and Plot #
params = linear_fit(thin1,hi1)
a      = params['a']
b      = params['b']
ea     = params['ea']
eb     = params['eb']


plt.subplot(2,1,1)
plt.plot(thin1, hi1, 'b^', ms=10, label='$N_{HI} vs N^{thin}_{HI}$')
plt.plot(thin1, a*np.array(thin1) + b, 'r-', linewidth=4, label='Best linear fit')

a  = round(a, 2)
b  = round(b, 2)
ea = round(ea, 2)
eb = round(eb, 2)

plt.grid()
plt.title('$N_{HI}$ vs $N^{thin}_{HI}$, SPONGE', fontsize = 35)
plt.ylabel('$N_{HI}$', fontsize = 35)
# plt.xlabel('$N^{thin}_{HI}$', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(10., 5., '$N_{HI}[10^{20} cm^{-2}] = ('+str(a)+'\pm'+str(ea) +')\cdot N^*_{HI}[10^{20} cm^{-2}] + ('+str(b)+'\pm'+str(eb)+')$', color='blue', fontsize=20)
plt.legend(loc='lower right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')

# Fit and Plot #
params = linear_fit(thin2,hi2)
a      = params['a']
b      = params['b']
ea     = params['ea']
eb     = params['eb']

plt.subplot(2,1,2)
plt.plot(thin2, hi2, 'b^', ms=10, label='$N_{HI} vs N^{thin}_{HI}$')
plt.plot(thin2, a*np.array(thin2) + b, 'r-', linewidth=4, label='Best linear fit')

a  = round(a, 2)
b  = round(b, 2)
ea = round(ea, 2)
eb = round(eb, 2)

plt.grid()
plt.title('$N_{HI}$ vs $N^{thin}_{HI}$, MS', fontsize = 35)
plt.ylabel('$N_{HI}$', fontsize = 35)
plt.xlabel('$N^{thin}_{HI}$', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(10., 5., '$N_{HI}[10^{20} cm^{-2}] = ('+str(a)+'\pm'+str(ea) +')\cdot N^*_{HI}[10^{20} cm^{-2}] + ('+str(b)+'\pm'+str(eb)+')$', color='blue', fontsize=20)
plt.legend(loc='lower right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')


plt.show()





## Fit and Plot #
## MS 79 sources
thin2 = np.asarray(thin2)
hi2   = np.asarray(hi2)
params = linear_fit(np.log10(thin2),hi2/thin2 )
a      = params['a']
b      = params['b']
ea     = params['ea']
eb     = params['eb']

plt.subplot(2,1,1)
plt.plot(np.log10(thin2),hi2/thin2, 'b^', ms=10, label='$N_{HI} vs N^{thin}_{HI}$')
plt.plot(np.log10(thin2), a*np.array(np.log10(thin2)) + b, 'r-', linewidth=4, label='Best linear fit')

a  = round(a, 2)
b  = round(b, 2)
ea = round(ea, 2)
eb = round(eb, 2)

plt.grid()
plt.title('$N_{HI}$ vs $N^{thin}_{HI}$, MS', fontsize = 35)
plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
# plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.2, 1.5, '$f = ['+str(a)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)
plt.legend(loc='lower right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')

## 21-SPONGE 30 src
thin1  = np.asarray(thin1)
hi1    = np.asarray(hi1)
params = linear_fit(np.log10(thin1),hi1/thin1 )
a      = params['a']
b      = params['b']
ea     = params['ea']
eb     = params['eb']

plt.subplot(2,1,2)
plt.plot(np.log10(thin1),hi1/thin1, 'b^', ms=10, label='$N_{HI} vs N^{thin}_{HI}$')
plt.plot(np.log10(thin1), a*np.array(np.log10(thin1)) + b, 'r-', linewidth=4, label='Best linear fit')

a  = round(a, 2)
b  = round(b, 2)
ea = round(ea, 2)
eb = round(eb, 2)

plt.grid()
plt.title('$N_{HI}$ vs $N^{thin}_{HI}$, 21-SPONGE', fontsize = 35)
plt.ylabel('$Ratio f = N_{HI}$/$N^*_{HI}$', fontsize=35)
plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.2, 1.5, '$f = ['+str(a)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='blue', fontsize=20)
plt.legend(loc='lower right', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')

plt.show()