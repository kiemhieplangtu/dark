import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl

import operator
from restore             import restore

## Read NHI from 94src #
 #
 # params string fname Filename
 # return dict info of N(HI)
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_nhi_94src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt'):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
	data = restore(fname, 2, cols, fmt)
	return data.read()


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

##================= MAIN ========================##
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
quantity = 'tau353'
map_unit = ''

# map_file = 'data/HFI_SkyMap_353_2048_R2.02_full.fits'
# quantity = 'w353'
# map_unit = 'K_CMB'


# info    = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
info  = read_nhi_94src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')

# Define the width of area #
beam   = 5.0            # Beam = 5' = map_resolution
dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
offset = dbeam          # degree

# Define constants #
nside     = 2048
deg2rad   = np.pi/180.

tau_map   = hp.read_map(map_file, field = 0)
nside     = hp.get_nside(tau_map)
res       = hp.nside2resol(nside, arcmin=False)
dd        = res/deg2rad/5.7

src = info['src']
for i in range(len(src)):
	if( (src[i] != '3C409') ):
		continue

	l  = info['l'][i]
	b  = info['b'][i]

	# Plot cartview a/o mollview #
	ll = l
	if (l>180):
		ll = ll-360.

	hp.cartview(tau_map, title=src[i], coord='G', unit='',
			norm=None, xsize=800, lonra=[ll-offset-0.1*offset,ll+offset+0.1*offset], latra=[b-offset-0.1*offset,b+offset+0.1*offset],
			return_projected_map=True)

	hp.projplot(ll, b, 'bo', lonlat=True, coord='G')
	hp.projtext(ll, b, ' (' + str(round(ll,2)) + ',' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=18, weight='bold', color='b')

	# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
	# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

	# Cal. #
	theta = (90.0-b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)

	for x in pl.frange(l-offset, l+offset, dd):
		for y in pl.frange(b-offset, b+offset, dd):
			cosb = np.cos(b*deg2rad)
			cosy = np.cos(y*deg2rad)
			if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
				# hp.projtext(x, y, '.', lonlat=True, coord='G')
				hp.projplot(x, y, 'kx', lonlat=True, coord='G', ms=12, mew=2)

	mpl.rcParams.update({'font.size':30})
	plt.show()