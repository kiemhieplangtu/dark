import os, sys
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import operator
import matplotlib.cm     as cm

from restore             import restore
from masks               import masks

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict infocd 
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_info_no_co(fname = '../gas_col_density/26src_no_co_info.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j','nhi','wnm','cnm','e_nhi','nhi_er','oh']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',   'f',  'f',  'f',  'f',    'f',    'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

#================= MAIN ========================#

## Classes ##
masks = masks()

## Define constants #
deg2rad = np.pi/180.
pi      = 2.0*np.arccos(0.)

## Read infor of 26 no-CO sources
info    = read_info_no_co('../../gas_col_density/26src_no_co_info.dat')

## All Masks
pth      = os.getenv("HOME")+'/hdata/dust/'
msk_file = pth + 'planck_mask.fits' #
msk      = hp.read_map(msk_file, field = 0, h=False)
title    = msk_file

nside    = hp.get_nside(msk)
# resol    = hp.nside2resol(nside, arcmin=False)
# print nside, resol*180./pi

cmap = plt.cm.get_cmap('spring')
hp.mollview(msk, title="Planck Mask", coord='G', unit='', rot=[0,0,0], norm=None, xsize=800, cmap=cmap)

for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	op,er = masks.get_dust_opacity(msk,l,b)
	hp.projplot(l, b, 'kx', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	hp.projtext(l, b, src+', '+str(op)+', '+str(er), lonlat=True, coord='G')

hp.graticule()
plt.show()