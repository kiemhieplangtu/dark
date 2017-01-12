import sys, os
sys.path.insert(0, r'/home/vnguyen/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
from   restore           import restore

## Read NHI from 94src #
 #
 # params string fname Filename
 # return dict info of N(HI)
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_nhi_94src(fname = '../result/nhi_thin_cnm_wnm_94src_sponge_prior.txt'):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
	data = restore(fname, 2, cols, fmt)
	return data.read()

#================= MAIN ========================#	
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/dust/LAB_Healpix/'
pth      = os.getenv("HOME")+'/hdata/hi/'
map_file = pth + 'LAB_fullvel.fits'

info     = read_nhi_94src(fname = '../result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
dbeam    = 30./120.0 # Beam = 0.6deg ~ 30' -> dbeam = beam/60/2
#dbeam = 0.1

## E(B-V) map ##
hi_map = hp.read_map(map_file, verbose=False, field = 0)
nside  = hp.get_nside(hi_map)
res    = hp.nside2resol(nside, arcmin = False)
dd     = res/deg2rad/2.0

#============= For Mask ========= #
offset = 2.0 #degree

#====== For Plotting ======#
hp.mollview(hi_map, title=r'$N_{HI}$', coord='G', unit='$cm^{-2}$', norm='hist') #, min=-1.6,max=1.1)

for i in range(len(info['src'])):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = hi_map[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'bo', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	if (b<60):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, weight='bold')
	else:
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='r', weight='bold')

	if(src == '3C109'):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='b', weight='bold')

mpl.rcParams.update({'font.size':30})
hp.graticule()
plt.grid()
plt.show()

sys.exit()

l  = 51.641648
b  = -9.6750019

# Plot cartview a/o mollview #
ll = l
if (l>180):
	ll = ll-360.

offset = 1.
lonr = [51.25, 52.0]
latr = [-10., -9.35]
hp.cartview(ci_map, title='XXX', coord='G', unit='mag', min=0.1,max=0.4,
		norm=None, xsize=800, lonra=lonr, latra=latr, #lonra=[ll-offset,ll+offset], latra=[b-offset,b+offset], min=0, max=0.4,
		return_projected_map=True)

# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

# Cal. #
hp.projplot(ll, b, 'bo', lonlat=True, coord='G')
hp.projtext(ll, b, ' (' + str(round(ll,2)) + ',' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=18, weight='bold')

# theta = (90.0-b)*deg2rad
# phi   = l*deg2rad
# pix   = hp.ang2pix(nside, theta, phi, nest=False)

# for x in pl.frange(l-offset, l+offset, dd):
# 	for y in pl.frange(b-offset, b+offset, dd):
# 		cosb = np.cos(b*deg2rad)
# 		cosy = np.cos(y*deg2rad)
# 		if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
# 			# hp.projtext(x, y, '.', lonlat=True, coord='G')
# 			hp.projplot(x, y, 'kx', lonlat=True, coord='G')

mpl.rcParams.update({'font.size':30})
plt.show()