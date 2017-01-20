import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import operator
import matplotlib.cm     as cm

from restore             import restore

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

## Get the Dust opacity and it Standard Deviation for each pixel (l,b) #
 #
 # params float l Gal-longitude
 # params float b Gal-latitude
 # params array msk Map of masks
 #
 # return list [f,er] Dust opacity and it Standard Deviation 
 # 
 # version 11/2016
 # Author Van Hiep ##
def get_dust_opacity(msk,l,b):	
	deg2rad = np.pi/180.
	nside   = hp.get_nside(msk)
	opac    = [8.4e-27, 7.1e-27, 6.8e-27, 6.5e-27 ]
	err     = [3.0e-27, 1.9e-27, 1.8e-27, 1.8e-27 ]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = msk[pix]
	return [ opac[int(val)], err[int(val)] ]

#================= MAIN ========================#
## Define constants #
deg2rad = np.pi/180.
pi      = 2.0*np.arccos(0.)

## Read infor of 26 no-CO sources
info    = read_info_no_co('26src_no_co_info.dat')

## g35+g45+g56 Mask
map_file = 'data/mask_g35g45g56.fits' #
msk      = hp.read_map(map_file, field = 0, h=False)
title    = map_file

nside    = hp.get_nside(msk)
# resol    = hp.nside2resol(nside, arcmin=False)
# print nside, resol*180./pi

cmap = plt.cm.get_cmap('cool')
hp.mollview(msk, title=map_file, coord='G', unit='', rot=[0,0,0], norm=None, xsize=800, cmap=cmap)

for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	op,er = get_dust_opacity(msk, l, b)
	hp.projplot(l, b, 'kx', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	hp.projtext(l, b, src+', '+str(op)+', '+str(er), lonlat=True, coord='G')

hp.graticule()
plt.show()