import sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

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
def read_nhi(fname = ''):
	cols = ['indx','src', 'nhi', 'nhi_er']
	fmt  = ['i',    's',   'f',  'f']
	data = restore(fname, 2, cols, fmt)
	return data.read()

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
info   = read_nhi('../result/total_nhi_1sig_uncertainty.txt')
srcs   = info['src']
nhi    = info['nhi']
nhi_er = info['nhi_er']
for i in range(len(srcs)):
	print ('{}    {}\t{}     {}     {}     {}'\
		.format(i, srcs[i],lbx[srcs[i]]['l'],lbx[srcs[i]]['b'], nhi[i], nhi_er[i] ))

plt.plot(b,nhi,'r.')
# plt.plot(l,nhi,'r.')
plt.grid()
plt.show()