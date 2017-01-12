import sys, os
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
from   restore           import restore

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
ms78sc = read_info_ms_78src(fname = '../result/nhi_lb_thin_cnm_wnm_78src.txt')
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

spnhi = []
msnhi = []
spwnm = []
mswnm = []
spcnm = []
mscnm = []
xl    = []
xb    = []
for i in range(len(comsc)):
	sc = comsc[i]
	# print sc, spdat[sc]['l'], msdat[sc]['l'], spdat[sc]['l'] - msdat[sc]['l']
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(spdat[sc]['nhi'] - msdat[sc]['nhi'])/msdat[sc]['nhi']
	spnhi.append(spdat[sc]['nhi'])
	msnhi.append(msdat[sc]['nhi'])
	spwnm.append(spdat[sc]['wnm'])
	mswnm.append(msdat[sc]['wnm'])
	spcnm.append(spdat[sc]['cnm'])
	mscnm.append(msdat[sc]['cnm'])
	xl.append(msdat[sc]['l'])
	xb.append(msdat[sc]['b'])

print 'NHI'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(-spdat[sc]['nhi'] + msdat[sc]['nhi'])/msdat[sc]['nhi']

print 'WNM'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['wnm'], msdat[sc]['wnm'], 100.*(-spdat[sc]['wnm'] + msdat[sc]['wnm'])/msdat[sc]['wnm']

print 'CNM'
for i in range(len(comsc)):
	sc = comsc[i]
	if(msdat[sc]['cnm'] != 0.):
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm'], 100.*(-spdat[sc]['cnm'] + msdat[sc]['cnm'])/msdat[sc]['cnm']
	else:
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm']

spnhi = np.asarray(spnhi)
msnhi = np.asarray(msnhi)
y     = 100.*(msnhi-spnhi)/msnhi
plt.plot(y, 'r*-', ms=12, label='$[(N^{MS}_{HI}-N^{SP}_{HI})/N^{MS}_{HI}] (\%)$')
plt.grid()
plt.title('$N_{HI}$ Difference, SPONGE vs MS', fontsize = 35)
plt.ylabel('$[(N^{MS}_{HI}-N^{SP}_{HI})/N^{MS}_{HI}] (\%)$', fontsize = 35)
plt.xlabel('Index', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.1, 0.,'14 common sources',color='b',fontsize=18)
plt.legend(loc='upper left', fontsize = 18)
# plt.axvline(x=60., lw=4)
for i in range(len(comsc)):
	sl = str(round(xl[i],2) )
	sb = str(round(xb[i],2) )
	plt.annotate('('+str(comsc[i])+') '+sl+', '+sb, xy=(i, y[i]), xycoords='data',
           xytext=(-50.,30.), textcoords='offset points',
           arrowprops=dict(arrowstyle="->"),fontsize=18,
           )
plt.show()

## WNM
plt.subplot(2,1,1)
spwnm = np.asarray(spwnm)
mswnm = np.asarray(mswnm)
y     = 100.*(mswnm-spwnm)/mswnm
plt.plot(y, 'r*-', ms=12, label='$[(N^{MS}_{wnm}-N^{SP}_{wnm})/N^{MS}_{wnm}] (\%)$')
plt.grid()
plt.title('$N_{wnm}$ Difference, SPONGE vs MS', fontsize = 35)
plt.ylabel('$[(N^{MS}_{wnm}-N^{SP}_{wnm})/N^{MS}_{wnm}] (\%)$', fontsize = 35)
# plt.xlabel('Index', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.1, 0.,'14 common sources',color='b',fontsize=18)
plt.legend(loc='upper left', fontsize = 18)
# plt.axvline(x=60., lw=4)
for i in range(len(comsc)):
	sl = str(round(xl[i],2) )
	sb = str(round(xb[i],2) )
	plt.annotate('('+str(comsc[i])+') '+sl+', '+sb, xy=(i, y[i]), xycoords='data',
           xytext=(-50.,-20.), textcoords='offset points',
           arrowprops=dict(arrowstyle="->"),fontsize=12,
           )

## CNM
plt.subplot(2,1,2)
spcnm = np.asarray(spcnm)
mscnm = np.asarray(mscnm)
y     = 100.*(mscnm-spcnm)/mscnm
plt.plot(y, 'r*-', ms=12, label='$[(N^{MS}_{cnm}-N^{SP}_{cnm})/N^{MS}_{cnm}] (\%)$')
plt.grid()
plt.title('$N_{cnm}$ Difference, SPONGE vs MS', fontsize = 35)
plt.ylabel('$[(N^{MS}_{cnm}-N^{SP}_{cnm})/N^{MS}_{cnm}] (\%)$', fontsize = 35)
plt.xlabel('Index', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.1, 0.,'14 common sources',color='b',fontsize=18)
plt.legend(loc='upper left', fontsize = 18)
# plt.axvline(x=60., lw=4)
for i in range(len(comsc)):
	sl = str(round(xl[i],2) )
	sb = str(round(xb[i],2) )
	plt.annotate('('+str(comsc[i])+') '+sl+', '+sb, xy=(i, y[i]), xycoords='data',
           xytext=(-50.,-20.), textcoords='offset points',
           arrowprops=dict(arrowstyle="->"),fontsize=12,
           )
plt.show()