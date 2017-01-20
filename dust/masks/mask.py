import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
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

## G35, G45, G56 masks, NSIDE = 2048 ##
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth+'COM_Mask_Likelihood_2048_R1.10.fits'
# TTYPE1  = 'CL31    '           /                                                
# TTYPE2  = 'CL39    '           /                                                
# TTYPE3  = 'CL49    '           /                                                
# TTYPE4  = 'G22     '           /                                                
# TTYPE5  = 'G35     '           /                                                
# TTYPE6  = 'G45     '           /                                                
# TTYPE7  = 'G56     '           /                                                
# TTYPE8  = 'G65     '           /                                                
# TTYPE9  = 'PS96    '           /                                                
# TTYPE10 = 'PSA82   '           /

## Other masks, NSIDE = 256 ##
lowest  = hp.read_map(pth+'mask_Lowest1perc_ns256.fits', field = 0, h=False)
low_hi  = hp.read_map(pth+'mask_LowNHI_ns256.fits',      field = 0, h=False)
sth_cap = hp.read_map(pth+'mask_SouthCap_ns256.fits',    field = 0, h=False)
gt15    = hp.read_map(pth+'mask_GLATgt15_ns256.fits',    field = 0, h=False)
whole   = hp.read_map(pth+'mask_WholeSky_ns256.fits',    field = 0, h=False)

## All Masks ##
lowest  = hp.ud_grade(lowest, nside_out=2048)    # weight = 7
low_hi  = hp.ud_grade(low_hi, nside_out=2048)    # weight = 6
sth_cap = hp.ud_grade(sth_cap,nside_out=2048)    # weight = 5
g35     = hp.read_map(map_file, field = 4, h=0)  # weight = 4
g45     = hp.read_map(map_file, field = 5, h=0)  # weight = 3
gt15    = hp.ud_grade(gt15,     nside_out=2048)  # weight = 2
g56     = hp.read_map(map_file, field = 6, h=0)  # weight = 1
whole   = hp.ud_grade(whole,    nside_out=2048)  # weight = 0

# msk   = whole
# msk   = g56
# msk   = gt15
# msk   = g45
# msk   = g35
# msk   = sth_cap
# msk   = low_hi
# msk   = lowest

msk = lowest+low_hi+sth_cap+g35+g45+gt15+g56

## Color map ##
nside = 2048
cmap  = plt.cm.get_cmap('spring')
cmap  = mpl.colors.ListedColormap(['c', 'r', 'g', 'b', 'y', 'm', 'k'])
hp.mollview(msk, title='Mask', coord='G', unit='', rot=[0,0,0], norm=None, xsize=800, cmap=cmap)

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
# hp.write_map("planck_mask.fits", msk)
plt.show()