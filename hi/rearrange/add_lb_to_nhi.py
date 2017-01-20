import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

## Read 78 sources with l,b and ra, dec
 # Note: 3C223 with No N(HI) and N(HI)_uncertainty #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 04/10/2016
 # Author Van Hiep ##
def read_lb(fname = ''):
	cols = ['l', 'b', 'src', 'ra', 'dec']
	fmt  = ['f', 'f',  's',   'f',  'f']
	data = restore(fname, 3, cols, fmt)
	return data.read()

## Read 78 sources with N(HI) and uncertainties
 # Note: 3C223 with No N(HI) and N(HI)_uncertainty #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 04/10/2016
 # Author Van Hiep ##
def read_nhi(fname = '../result/78src_nhi_cnm_wnm_with_err.txt'):
	cols = ['idx','src',  'nhi', 'nhi_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',    's',   'f',  'f',       'f',   'f',       'f',   'f' ]
	data = restore(fname, 2, cols, fmt)
	return data.read()

## Read 78 sources with N(HI)_thin and its uncertainties
 # Note: 3C223 with No N(HI) and N(HI)_uncertainty #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_nhi_thin(fname = '../result/nhi_thin_with_er_78src.txt'):
	cols = ['idx','src', 'l', 'b', 'thin', 'thin_er']
	fmt  = ['i',    's',  'f', 'f', 'f',  'f']
	data = restore(fname, 2, cols, fmt)

	info = data.read()
	src  = info['src']
	thin = info['thin']
	err  = info['thin_er']

	ret = {}
	for i in range(len(src)):
		sc              = src[i] 

		ret[sc]         = {}
		ret[sc]['thin'] = thin[i]
		ret[sc]['er']   = err[i]

	return ret

##=================== MAIN ===================##
lb  = read_lb('../data/78src_radec_lb.txt')
src = lb['src']
l   = lb['l']
b   = lb['b']

lbx = {}
for i in range(len(src)):
	lbx[src[i]] = {}
	lbx[src[i]]['l'] = lb['l'][i]
	lbx[src[i]]['b'] = lb['b'][i]

## Read N(HI) ##
thin = read_nhi_thin(fname = '../result/nhi_thin_with_er_78src.txt')

## Read N(HI) ##
info   = read_nhi('../result/78src_nhi_cnm_wnm_with_err.txt')
srcs   = info['src']
nhi    = info['nhi']
nhi_er = info['nhi_er']
cnm    = info['cnm']
cnm_er = info['cnm_er']
wnm    = info['wnm']
wnm_er = info['wnm_er']
for i in range(len(srcs)):
	sc = srcs[i]
	print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
		# .format(i, srcs[i],lbx[srcs[i]]['l'],lbx[srcs[i]]['b'], nhi[i], nhi_er[i] ))                   ## file: hi/rearrange/nhi_lb_78src.txt
		.format(i, sc,lbx[sc]['l'],lbx[sc]['b'], nhi[i], nhi_er[i], thin[sc]['thin'], thin[sc]['er'], cnm[i], cnm_er[i], wnm[i], wnm_er[i]  ))  ## file: hi/rearrange/nhi_lb_thin_cnm_wnm_78src.txt

plt.plot(b,nhi,'r.')
# plt.plot(l,nhi,'r.')
plt.grid()
plt.show()