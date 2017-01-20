import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp

import operator
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

#================= MAIN ========================#	
## Define constants ##
deg2rad = np.pi/180.
beam    = 3.5
dbeam   = beam/120.0 # Beam = 3.5' -> dbeam = beam/60/2
dbeam   = 0.5
edge    = dbeam/40. # To plot a rectangle about the source
info    = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')

## CO Dame Map, NSIDE = 512 ##
map_file = os.getenv("HOME")+'/hdata/co/lambda_wco_dht2001.fits'
co_map   = hp.read_map(map_file, field = 1, h=False)

## Color map ##
nside = 512
# cmap  = plt.cm.get_cmap('winter')
hp.mollview(co_map, title='Velocity Integrated CO Map', coord='G', unit='K km/sec', rot=[0,0,0], norm='hist', xsize=800)

for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = msk[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'bx', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	hp.projtext(l, b, src+','+str(val), lonlat=True, coord='G')

hp.graticule()
plt.show()