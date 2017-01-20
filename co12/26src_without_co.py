import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import pylab             as pl
import operator

from numpy    import array
from restore  import restore

## Read Basic info of 26 sources #
 #
 # params string fname Filename
 # return void
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_info_no_co(fname = 'result/26src_no_co_basic.dat'):
	cols = ['id','src','l','b','ra_icrs','de_icrs','ra_j','de_j','oh']
	fmt  = ['i', 's',  'f','f','f',       'f',     's',    's',  'i']
	data = restore(fname, 2, cols, fmt)
	return data.read()

## Read NHI from 94src #
 #
 # params string fname Filename
 # return dict info of N(HI)
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_nhi_94src(fname = '../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt'):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
	data = restore(fname, 2, cols, fmt)
	return data.read()

## ========= MAIN ============= ##
inf26 = read_info_no_co(fname = 'result/26src_no_co_basic.dat')
sc26  = inf26['src']
l     = inf26['l']
b     = inf26['b']
rai   = inf26['ra_icrs']
deci  = inf26['de_icrs']
raj   = inf26['ra_j']
decj  = inf26['de_j']
oh    = inf26['oh']

## 94 HI sources
inf94  = read_nhi_94src(fname = '../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
sc94   = inf94['src']
nhi    = inf94['nhi']
nhier  = inf94['nhi_er']
thin   = inf94['thin']
thiner = inf94['thin_er']
cnm    = inf94['cnm']
cnmer  = inf94['cnm_er']
wnm    = inf94['wnm']
wnmer  = inf94['wnm_er']

for i in range(len(sc26)):
	k = sc94.index(sc26[i])
	print('{:3}  {:11} {:.4f}    {:.4f}    {:.4f}    {:.4f}    {}    {}    {}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}'\
		.format(i, sc26[i], l[i], b[i], rai[i], deci[i], raj[i], decj[i], oh[i], nhi[k], nhier[k], thin[k], thiner[k], cnm[k], cnmer[k], wnm[k], wnmer[k] ) )      ## result/26src_no_co_with_sponge.dat