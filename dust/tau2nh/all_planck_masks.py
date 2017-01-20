import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

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
 # version 1/2017
 # Author Van Hiep ##
def read_info_no_co(fname = '../../co12/result/26src_no_co_with_sponge.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j', 'oh', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',    'i', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
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
info    = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')

## g35+g45+g56 Mask
msk_file = '../gas_col_density/data/mask_WholeSky_ns256.fits' #
msk_file = '../gas_col_density/data/mask_GLATgt15_ns256.fits' #
msk_file = '../gas_col_density/data/mask_Lowest1perc_ns256.fits' #
msk_file = '../gas_col_density/data/mask_LowNHI_ns256.fits' #
# msk_file = '../gas_col_density/data/mask_LVC3E20_ns256.fits' #
# msk_file = '../gas_col_density/data/mask_SouthCap_ns256.fits' #
msk      = hp.read_map(msk_file, field = 0, h=False)
msk      = hp.ud_grade(msk,nside_out=2048)
title    = msk_file

nside    = hp.get_nside(msk)
print 'Nside',nside
# resol    = hp.nside2resol(nside, arcmin=False)
# print nside, resol*180./pi

cmap = plt.cm.get_cmap('cool')
hp.mollview(msk, title=msk_file, coord='G', unit='', rot=[0,0,0], norm=None, xsize=800, cmap=cmap)


# for i in range(0,26):
# 	src   = info['src'][i]
# 	l     = info['l'][i]
# 	b     = info['b'][i]

# 	op,er = masks.get_dust_opacity(msk,l,b)
# 	hp.projplot(l, b, 'kx', lonlat=True, coord='G')
# 	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
# 	hp.projtext(l, b, src+', '+str(op)+', '+str(er), lonlat=True, coord='G')

hp.graticule()
plt.show()