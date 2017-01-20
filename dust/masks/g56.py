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

def discrete_cmap(N, base_cmap=None):
	base = plt.cm.get_cmap(base_cmap)
	color_list = base(np.linspace(0, 1, N))
	cmap_name = base.name + str(N)
	return base.from_list(cmap_name, color_list, N)

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

## G56 Mask ##
msk   = hp.read_map(map_file, field = 4, h=0)  # weight = 1

## Color map ##
nside = 2048
cmap  = plt.cm.get_cmap('spring')
cmap  = mpl.colors.ListedColormap(['lightgray', 'firebrick'])
hp.mollview(msk, title='G35', coord='G', unit='', rot=[0,0,0], norm=None, xsize=800, cmap=cmap)

hp.graticule()
# hp.write_map("planck_mask.fits", msk)
plt.show()